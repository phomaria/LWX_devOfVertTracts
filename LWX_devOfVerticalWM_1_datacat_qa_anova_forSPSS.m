% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';

% beh_measure = 'age'; %age, lit, vm, fm
wm_measure_here = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'}; 

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure_here{w} '_raw.mat']))
    
    % Convert into array and header for ease.
    data_all_in = table2array(data_tbl);
    data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Get index matrices for hypothesis-driven grouping of WM tracts. 
    for k = 1:length(data_all_in_header)
        
        % Indices of horizontal tracts.
        h_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
            || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF') ...
            || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
            || strcmp(data_all_in_header{k}, 'leftArc') || strcmp(data_all_in_header{k}, 'rightArc') ...
            || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3');
        
        % Indices of vertical tracts.
        v_idx(k) = strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
            || strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
            || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
            || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang');
        
        % Set the grouping variable for horizontal (=1) and vertical (=2) tracts and tracts that are not of interest (=0).
        if h_idx(k) == 1
            
            hv(k) = 1;
            
        elseif v_idx(k) == 1
            
            hv(k) = 2;
            
        else
            
            hv(k) = 0;
            
        end
        
    end
    
    % Select the measurements of the tracts that I care about and convert all zeros to NaN.
    toi = data_all_in(:, find(hv ~= 0)); toi(toi==0) = NaN;
    
    % Update the tract list so that we can grab the tract name.
    data_all_in_header = data_all_in_header(hv~=0);
    
    % Update the tract indexing to correspond with toi dimensions.
    hv = hv(find(hv~=0));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(toi(find(data_tbl.gp_age == 1), :));
    m_wtwg2 = nanmean(toi(find(data_tbl.gp_age == 2), :));
    m_wtwg3 = nanmean(toi(find(data_tbl.gp_age == 3), :));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(toi(find(data_tbl.gp_age == 1), :));
    std_wtwg2 = nanstd(toi(find(data_tbl.gp_age == 2), :));
    std_wtwg3 = nanstd(toi(find(data_tbl.gp_age == 3), :));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2;
    r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3;
    r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, [size(toi(find(data_tbl.gp_age==1)))]), repmat(r_max_wtwg2, size(toi(find(data_tbl.gp_age==2)))), ...
        repmat(r_max_wtwg3, size(toi(find(data_tbl.gp_age==3)))));
    r_min = cat(1, repmat(r_min_wtwg1, [size(toi(find(data_tbl.gp_age==1)))]), repmat(r_min_wtwg2, size(toi(find(data_tbl.gp_age==2)))), ...
        repmat(r_min_wtwg3, size(toi(find(data_tbl.gp_age==3)))));
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(toi > r_max))
        toi(find(toi > r_max)) = NaN;
        disp(['Replaced ' num2str(numel(find(toi > r_max))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.'])
    else
        disp('No data points were above 3 standard deviations of the within-tract, within-group mean.')
    end
    %Min
    if ~isempty(find(toi < r_min))
        toi(find(toi > r_min)) = NaN;
        disp(['Replaced ' num2str(numel(find(toi > r_min))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.'])
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.')
    end
    
    % Output csv file for ANOVA in SPSS. (Matlab doesn't handle Mixed Model
    % ANOVAs well when the between-group variable is correlated with subID
    % (e.g., when between-group variable is something like age groups).
    t_out = array2table(cat(2, data_tbl.subID, data_tbl.gp_age, data_tbl.gp_lit, data_tbl.gp_vm, ...
        data_tbl.gp_fm, data_tbl.cov_age, data_tbl.cov_sex, toi(:, 1:end), nanmean(toi(:, hv == 1), 2), nanmean(toi(:, hv == 2), 2)), 'VariableNames', ...
        {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', data_all_in_header{:}, 'meanH', 'meanV'});
    
    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '.csv']));
    
    % Output z-scored file for SPSS.
    toi_z = (nanmean(toi, 1) - toi)./nanstd(toi, [], 1);
    temp = nanmean(toi(:, hv == 1), 2);
    toi_meanh_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    temp = nanmean(toi(:, hv == 2), 2);
    toi_meanv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    t_out_z = array2table(cat(2, data_tbl.subID, data_tbl.gp_age, data_tbl.gp_lit, data_tbl.gp_vm, ...
        data_tbl.gp_fm, data_tbl.cov_age, data_tbl.cov_sex,toi_z(:, 1:end), toi_meanh_z, toi_meanv_z), 'VariableNames', ...
        {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', data_all_in_header{:}, 'meanH', 'meanV'});
    writetable(t_out_z, fullfile(rootDir, 'supportFiles', ['LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '_z.csv']));
    
end


