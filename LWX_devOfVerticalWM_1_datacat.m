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

list_wm_toi = {'leftAslant', 'rightAslant', 'leftVOF', 'rightVOF', 'leftTPC', 'rightTPC', ...
    'leftMDLFspl', 'rightMDLFspl', 'leftMDLFang', 'rightMDLFang', 'leftpArc', 'rightpArc', ...
    'leftArc', 'rightArc', 'leftSLF1And2', 'rightSLF1And2', 'leftCST', 'rightCST'};

remove_outliers = 'yes';

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
addpath(genpath([rootDir blprojectid]));

% Read in behavioral data.
load([rootDir 'supportFiles/LWX_all_groupings.mat']);

% Header for LWX Behavioral DatFa. gp assignments are
% young children (=1), older children (=2), and adults (=3).
% gp_lit (early literate = 1, literate = 2),
% gp_vm (low cm = 1, high vm = 2), and
% gp_fm (low cm = 1, high vm = 2).
header = {'subID', 'age_mo', 'vmi', 'vp', 'mc', 'pegs_dom', 'pegs_ndom', 'lwi', 'spell', 'wa', 'sos', 'c_vm', 'c_fm', 'c_lit', 'gp_age', 'gp_lit', 'gp_vm', 'gp_fm'};

fcount = 0;
for w = 1%:length(w_measures)
    
    for b = 1%:length(b_measures)
        
        wm_measure = w_measures{w};
        beh_measure = b_measures{b};
        
        %% BEHAVIOR.
        
        % Set names for correct columns given user input..
        if strcmp(beh_measure, 'age')
            
            measure_in = 'age_mo';
            gp_in = 'gp_age';
            
            outlier = [128 318];
            %128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
            % 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
            
        elseif strcmp(beh_measure, 'lit')
            
            measure_in = 'c_lit';
            gp_in = 'gp_lit';
            
            outlier = [128 318];
            
        elseif strcmp(beh_measure, 'vm')
            
            measure_in = 'c_vm';
            gp_in = 'gp_vm';
            
            outlier = [128 318];
            
        elseif strcmp(beh_measure, 'fm')
            
            measure_in = 'c_fm';
            gp_in = 'gp_fm';
            
            outlier = [128 318];
            
        end
        
        %% TRACTOGRAPHY.
        
        % Get contents of the directory where the tract measures for this subject are stored.
        grp_contents = dir([rootDir filesep blprojectid filesep]);
        
        % Remove the '.' and '..' files.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
        
        % Keep only names that are subject folders.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
        
        % Load in each tract's tractography measures for this subject.
        sub_count = 0; tract = {};
        for i = 1:size(grp_contents, 1)
            
            %     % Only read in data for subjects that we want to consider in this
            %     analysis. Right now we don't have all behavioral data entered.
            if sum(ismember(data_lwx(:, 1), grp_contents(i).name)) ~= 0
                
                % Display current sub ID.
                disp(grp_contents(i).name)
                
                % Update subject counter for when not all subjects are used/needed.
                sub_count = sub_count + 1;
                
                % Get contents of the directory where the tract measures for this subject are stored.
                sub_contents_tractprofiles = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-tractprofile*/profiles/*.csv']);
                sub_contents_tractstats = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-tractmeasures*/*.csv']);
                
                % Get contents of the directory where the SNR values for this subject are stored.
                sub_contents_snr = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-raw.tag-snr*/*product.json']);
                
                % Remove the '.' and '..' files.
                sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
                sub_contents_tractstats = sub_contents_tractstats(arrayfun(@(x) x.name(1), sub_contents_tractstats) ~= '.');
                sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
                
                % Read in stats for this subject's tracts.
                data_stats_temp = readtable([sub_contents_tractstats.folder filesep sub_contents_tractstats.name]);
                
                for j = 1:size(sub_contents_tractprofiles)
                    
                    % Read in data for this subject and this tract.
                    data_temp = csvread([sub_contents_tractprofiles(j).folder filesep sub_contents_tractprofiles(j).name], 1, 0);
                    
                    % Get middle 80%.
                    start = size(data_temp, 1)*.1;
                    stop = size(data_temp, 1)*.9;
                    
                    % Read in mean WM measure.
                    if strcmp(wm_measure, 'ad')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 1)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 2));
                        
                    elseif strcmp(wm_measure, 'fa')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 3)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 4));
                        
                    elseif strcmp(wm_measure, 'md')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 5)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 6));
                        
                    elseif strcmp(wm_measure, 'rd')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 7)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 8));
                        
                    elseif strcmp(wm_measure, 'icvf')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 17)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 18));
                        
                    elseif strcmp(wm_measure, 'isovf')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 19)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 20));
                        
                    elseif strcmp(wm_measure, 'od')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 21)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 22));
                        
                    end
                    
                    % Grab tract name for grouping variable.
                
                    tract{sub_count, j} = sub_contents_tractprofiles(j).name(1:end-13)
                    
                    % Grab subID.
                    sub(sub_count) = str2num(grp_contents(i).name(end-2:end));
                    
                    clear data_temp
                    
                    % Deal with differences in puncuation, spelling, etc.
                    if strcmp(tract{sub_count, j}, 'left_AnterioFrontoCerebellar')
                        
                        name_temp = strcat(tract{sub_count, j}(1:4), [' ' tract{sub_count, j}(6:end)]);
                        
                    elseif strcmp(tract{sub_count, j}, 'right_AnterioFrontoCerebellar')
                        
                        name_temp = strcat(tract{sub_count, j}(1:5), [' ' tract{sub_count, j}(7:end)]);
                        
                    elseif strncmp(tract{sub_count, j}, 'left_', 5)
                        
                        name_temp = strcat(tract{sub_count, j}(1:4), tract{sub_count, j}(6:end));
                        
                    elseif strncmp(tract{sub_count, j}, 'right_', 6)
                        
                        name_temp = strcat(tract{sub_count, j}(1:5), tract{sub_count, j}(7:end));
                        
                    else
                        
                        name_temp = tract{sub_count, j};
                        
                    end
                    
                    % Get number of streamlines for this tract and for this subject.
                    n_streamlines(sub_count, j) = data_stats_temp.StreamlineCount(find(strcmp(data_stats_temp.TractName, name_temp)));
                    
                    clear name_temp
                    
                end % end j
                
            end % ismember
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            for g = 1:size(data_snr_temp.SNRInB0_X_Y_Z, 1)
                
                get_temp(g) = str2num(data_snr_temp.SNRInB0_X_Y_Z{g});
                
            end
            
            snr(sub_count) = min(get_temp);
            
            clear data_stats_temp data_snr_temp get_temp
            
        end
        
        
    end
    
    % Set the mean of any tracts with less than 100 streamlines to NaN so that they will not be included in the average.
%     m(n_streamlines <= 100) = NaN;
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get WM measurements for each tract (reorganizing so that each column is a tract).
    for k = 1:size(list_tract, 1)
        
        % Get the mean of the wm_measure for this tract (take sum because don't want to include zero columns; only one value of interest per row).
        y_temp(:, k) = nansum(m.*strcmp(tract, list_tract{k}), 2);
        
    end % end k
    y = cat(2, transpose(sub), y_temp); clear y_temp
    
    % Create new grouping vectors, one subject at a time.
    for s = 1:length(sub)
        
        if ismember(sub(s), data_lwx(:, find(strcmp(header, 'subID'))))
            
            % Get this subject's AGE group.
            group(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, gp_in)));
            
            % Get this subject's MEASURE.
            measure(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, measure_in)));
            
            % Get this subject's AGE.
            cov_age(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, 'age_mo')));
            
            % Get this subject's WM measurements.
            wm(s, :) = y(find(y(:, 1) == sub(s)), 3:end); % start at 3 because first column is sub and second column is 'empty' (i.e., 0)
            
        end
        
    end % end s
    
    % Use only subjects for whom we have both behavioral and tractography measures.
    
    % Remove the 'group' of zero (i.e., empty). Make sure to modify 'group' last.
    measure = transpose(measure(group~=0));
    n_streamlines = n_streamlines(group~=0, :);
    wm = wm(group~=0, :);
    cov_age = transpose(cov_age(group~=0));
    sub = transpose(sub(group~=0));
    snr = transpose(snr(group~=0));
    group = transpose(group(group~=0));
    
    % Set all zeros to NaN.
    measure(measure == 0) = NaN; 
    wm(wm == 0) = NaN; 
    cov_age(cov_age == 0) = NaN; 
    sub(sub == 0) = NaN; 
    n_streamlines(n_streamlines == 0) = NaN; 
    snr(snr == 0) = NaN; 
    group(group == 0) = NaN;
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        idx_remove = ~any(~(sub~=outlier), 2);
        
        wm = wm(idx_remove, :);
        n_streamlines = n_streamlines(idx_remove, :);
        measure = measure(idx_remove);
        cov_age = cov_age(idx_remove);
        group = group(idx_remove);
        snr = snr(idx_remove);
        sub = sub(idx_remove);
        
    end
    
    % Make a data set for only the child data.
    wm_childrenOnly = wm(group~=3, :);
    n_streamlines_childrenOnly = n_streamlines(group~=3, :);
    measure_childrenOnly = measure(group~=3);
    cov_age_childrenOnly = cov_age(group~=3);
    sub_childrenOnly = sub(group~=3);
    snr_childrenOnly = snr(group~=3);
    group_childrenOnly = group(group~=3);
    
    % Save all variables.
    save([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_tractz.mat'])
    
    % Reset for next loop.
    clearvars -except w b rootDir header fcount blprojectid b_measures w_measures list_wm_toi test_type pval p_bonf remove_outliers save_figures data_lwx
    
end








