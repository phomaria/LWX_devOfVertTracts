% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

prefix = ''; % 'streamlinecount25920-'
blprojectid = ['/' prefix 'proj-5a74d1e26ed91402ce400cca/'];

b_measures = {'age', 'lit', 'vm', 'fm'};
w_measures = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'};

remove_outliers = 'yes';

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
    
    for b = 1:length(b_measures)
        
        wm_measure = w_measures{w};
        beh_measure = b_measures{b};       
        
        %% BEHAVIOR.
        
        % Set names for correct columns given user input..
        if strcmp(beh_measure, 'age')
            
            measure_in = 'Age_months';
            gp_in = 'group_age';
            
        elseif strcmp(beh_measure, 'lit')
            
            measure_in = 'c_lit';
            gp_in = 'group_lit';
                        
        elseif strcmp(beh_measure, 'vm')
            
            measure_in = 'c_vm';
            gp_in = 'group_vm';
                        
        elseif strcmp(beh_measure, 'fm')
            
            measure_in = 'c_fm';
            gp_in = 'group_fm';
                        
        end
        
        %% TRACTOGRAPHY.
        
        % Get contents of the directory where the tract measures for this subject are stored.
        grp_contents = dir([rootDir filesep blprojectid filesep]);
        
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
            sub_contents_tractprofiles = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-tractprofile*/profiles/*.csv']);
                        
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
                    
                    m(i, j) = mean(data_temp.ad_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'fa')
                    
                    m(i, j) = mean(data_temp.fa_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'md')
                    
                    m(i, j) = mean(data_temp.md_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'rd')
                    
                    m(i, j) = mean(data_temp.rd_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'icvf')
                    
                    m(i, j) = mean(data_temp.icvf_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'isovf')
                    
                    m(i, j) = mean(data_temp.isovf_1(start:stop)); 
                    
                elseif strcmp(wm_measure, 'od')
                    
                    m(i, j) = mean(data_temp.od_1(start:stop)); 
                    
                end
                
                % Grab tract name for grouping variable.
                tract{i, j} = sub_contents_tractprofiles(j).name(1:end-13);
                
                clear data_temp
                
            end % end j
            
        end
        
    end
    
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
    
    % Remove 'empty' column from data and header.
    wm = wm(:, find(~all(isnan(wm), 1)));
    wm_header = list_tract(~strcmp(list_tract, 'empty'));
    
    % Create grouping vectors, one subject at a time.
    for s = 1:length(sub)
        
        % Only do this for subjects for whom we have behavioral data.
        if ismember(sub(s), beh_data_in(:, find(strcmp(beh_data_in_header, 'SubjectID'))))
                    
            % Get this subject's AGE group.
            group(s) = beh_data_in(find(sub(s) == beh_data_in(:, find(strcmp(beh_data_in_header, 'SubjectID')))), find(strcmp(beh_data_in_header, gp_in)));
            
            % Get this subject's MEASURE.
            measure(s) = beh_data_in(find(sub(s) == beh_data_in(:, find(strcmp(beh_data_in_header, 'SubjectID')))), find(strcmp(beh_data_in_header, measure_in)));
            
            % Get this subject's AGE.
            cov_age(s) = beh_data_in(find(sub(s) == beh_data_in(:, find(strcmp(beh_data_in_header, 'SubjectID')))), find(strcmp(beh_data_in_header, 'Age_months')));
            
            % Get this subject's SEX.
            cov_sex(s) = beh_data_in(find(sub(s) == beh_data_in(:, find(strcmp(beh_data_in_header, 'SubjectID')))), find(strcmp(beh_data_in_header, 'Sex')));
                   
        end % end if ismember
        
    end % end s
          
    % Concatenate into one data array and one header array.
    data_all = cat(2, transpose(sub), transpose(group), transpose(measure), transpose(cov_age), transpose(cov_sex), wm); 
    data_all_header = [{'subID', gp_in, beh_measure, 'cov_age', 'cov_sex'}, wm_header{:}];
        
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(data_all(:, find(strcmp(data_all_header, {'subID'}))), outlier);
        
        % Remove outliers.
        data_all = data_all(~idx_outlier, :);
        
    end
        
    data_tbl = array2table(data_all, 'VariableNames', data_all_header);
    
    % Save all variables.
    save([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_raw.mat'], 'data_tbl')
    
    % Reset for next loop.
    clearvars -except w b rootDir beh_data_in_header beh_data_in blprojectid remove_outliers b_measures w_measures
    
end








