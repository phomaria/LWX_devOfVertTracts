% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

remove_outliers = 'yes';

prefix = ''; % 'streamlinecount25920-'
blprojectid = ['/' prefix 'proj-5a74d1e26ed91402ce400cca/'];

w_measures = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'};

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
addpath(genpath([rootDir blprojectid]));

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Parse table into array and header to make things easier.
beh_data_in = table2array(beh_data_in_tbl);
beh_data_in_header = beh_data_in_tbl.Properties.VariableNames;

% Identify outliers to be removed.
% 128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
% 315 because strange WM in occipital lobe leading to no right or left VOF
% 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
outlier = [128 315 318];

for w = 1:length(w_measures)
    
    wm_measure = w_measures{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    for i = 1:size(grp_contents, 1)
        
        % Grab subID.
        sub(i) = str2num(grp_contents(i).name(end-2:end));
        
        % Display current sub ID.
        disp(grp_contents(i).name)
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '/dt-neuro-tractprofile*/profiles/*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        for j = 1:size(sub_contents_tractprofiles)
            
            % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
            if i == 1 && j == 1
                
                tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                
            end
            
            % Read in data for this subject and this tract.
            data_temp = readtable([sub_contents_tractprofiles(j).folder filesep sub_contents_tractprofiles(j).name]);
            
            % Get middle 80%.
            start = size(data_temp, 1)*.1;
            stop = size(data_temp, 1)*.9;
            
            % Read in mean WM measure.
            if strcmp(wm_measure, 'ad')
                
                m(i, j) = nanmean(data_temp.ad_1(start:stop));
                
            elseif strcmp(wm_measure, 'fa')
                
                m(i, j) = nanmean(data_temp.fa_1(start:stop));
                
            elseif strcmp(wm_measure, 'md')
                
                m(i, j) = nanmean(data_temp.md_1(start:stop));
                
            elseif strcmp(wm_measure, 'rd')
                
                m(i, j) = nanmean(data_temp.rd_1(start:stop));
                
            elseif strcmp(wm_measure, 'icvf')
                
                m(i, j) = nanmean(data_temp.icvf_1(start:stop));
                
            elseif strcmp(wm_measure, 'isovf')
                
                m(i, j) = nanmean(data_temp.isovf_1(start:stop));
                
            elseif strcmp(wm_measure, 'od')
                
                m(i, j) = nanmean(data_temp.od_1(start:stop));
                
            end
            
            % Grab tract name for grouping variable.
            tract{i, j} = sub_contents_tractprofiles(j).name(1:end-13);
            
            clear data_temp
            
        end % end j
        
    end % end i
    
    % Find empty cells.
    t = find(cellfun(@isempty,tract));
    
    % Enter 'empty' in empty cells.
    tract(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get WM measurements for each tract (reorganizing so that each column is a tract).
    for k = 1:size(list_tract, 1)
        
        % Select the wm_measurements for this tract from each subject.
        temp = m.*strcmp(tract, list_tract{k});
        
        % Convert all zeros to NaN;
        temp(temp == 0) = NaN;
        
        % Get the mean of the wm_measure for this tract (take sum because don't want to include zero columns; only one value of interest per row).
        wm(:, k) = nansum(temp, 2);
        
        clear temp
        
    end % end k
    
    % Convert all zeros to NaN;
    wm(wm == 0) = NaN;
    
    % Remove 'empty' column from data and header and append subID.
    wm = cat(2, transpose(sub), wm(:, find(~all(isnan(wm), 1))));
    wm_header = [{'subID'}, transpose(list_tract(~strcmp(list_tract, 'empty')))];
    
    % Create grouping and behavioral vectors.
    beh = cat(2, beh_data_in_tbl.SubjectID, beh_data_in_tbl.group_age, beh_data_in_tbl.Age_months, beh_data_in_tbl.group_lit, beh_data_in_tbl.c_lit, ...
        beh_data_in_tbl.group_vm, beh_data_in_tbl.c_vm, beh_data_in_tbl.group_fm, beh_data_in_tbl.c_fm, beh_data_in_tbl.Sex);
    
    % Determine which subIDs appear in both WM and BEH.
    sub_wm_beh = intersect(wm(:, find(strcmp(wm_header, 'subID'))), beh_data_in_tbl.SubjectID);
        
    % Get indices of subjects who appear in both WM and BEH.
    sub_idx_wm = ismember(wm(:, find(strcmp(wm_header, 'subID'))), sub_wm_beh);
    sub_idx_beh = ismember(beh_data_in_tbl.SubjectID, sub_wm_beh);

    % Select only subjects who appear in both WM and BEH.
    % Concatenate into one data array and one header array.
    % Remove redundant subID columns.
    data_all = cat(2, beh(sub_idx_beh, :), wm(sub_idx_wm, find(strcmp(wm_header, 'subID'))+1:end));    
    data_all_header = [{'subID',  'gp_age', 'cov_age', 'gp_lit', 'c_lit', 'gp_vm', 'c_vm', ...
        'gp_fm', 'c_fm', 'cov_sex'}, wm_header{find(strcmp(wm_header, 'subID'))+1:end}];
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(data_all(:, find(strcmp(data_all_header, {'subID'}))), outlier);
        
        % Remove outliers.
        data_all = data_all(~idx_outlier, :);
        
    end
    
    data_tbl = array2table(data_all, 'VariableNames', data_all_header);
    
    % Save all variables.
    save([rootDir 'supportFiles/LWX_data_' wm_measure '_raw_singleshell.mat'], 'data_tbl')
    
    % Reset for next loop.
    clearvars -except w rootDir beh_data_in_tbl beh_data_in_header beh_data_in blprojectid remove_outliers w_measures outlier
    
end









