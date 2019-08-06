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

beh_measure = 'age'; %age, lit, vm, fm
wm_measure_here = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'}; %fa, ad, md, rd, od, icvf, isovf
sub = [103, 105, 108, 109, 115, 119, 121, 123, 125, 126, 127, 204, 213, 215, 216, ...
   219, 220, 221, 222, 226, 301, 302, 303, 304, 306, 314, 315, 316, 319];

% %% BEHAVIORAL MEASURES
% 
% % Read in behavioral data.
% load([rootDir 'supportFiles/LWX_all_groupings.mat']);
% 
% % Header for LWX Behavioral Data.
% header = {'subID', 'age_mo', 'vmi', 'vp', 'mc', 'pegs_dom', 'pegs_ndom', 'lwi', 'spell', 'wa', 'sos', 'c_vm', 'c_fm', 'c_lit', 'gp_age', 'gp_lit', 'gp_vm', 'gp_fm'};
% 
% % Select only those subjects for whom we have both behavioral and tractography data.
% subID_list = data_lwx(find(ismember(sub, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'subID')));
% 
% % Subselect only data for children and for the covariates of interest. Remove unnecessary columns.
% beh = cat(2, data_lwx(find(ismember(sub, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'age_mo'))), ...
%     data_lwx(find(ismember(sub, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'c_lit'))), ...
%     data_lwx(find(ismember(sub, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'c_vm'))), ...
%     data_lwx(find(ismember(sub, data_lwx(:, strcmp(header, 'subID')))), find(strcmp(header, 'c_fm'))));
% 
% % Get measure-specific z-scores (even for age).
% z_beh = (beh - nanmean(beh))./nanstd(beh);

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load([rootDir 'supportFiles/LWX_data_' wm_measure_here{w} '_' beh_measure '_tractz_test32104.mat'])
    
    % Tidy-up working directory.
    clearvars -except w rootDir beh_measure wm_measure_here list_tract sub cov_age group wm beh z_beh
    
    % Get index matrices for hypothesis-driven grouping of WM tracts. LH only right now.
    for k = 1:size(list_tract, 1)
        
        % Indices of horizontal tracts.
        h_idx(k) = strcmp(list_tract{k}, 'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
            || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
            || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
            || strcmp(list_tract{k}, 'leftArc') || strcmp(list_tract{k}, 'rightArc') ...
            || strcmp(list_tract{k}, 'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3');
        
        % Indices of vertical tracts.
        v_idx(k) = strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
            || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
            || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
            || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
            || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
            || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang');
        
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
    toi = wm(:, find(hv ~= 0)); toi(toi==0) = NaN;
    
    % Update the tract list so that we can grab the tract name.
    list_tract = list_tract(hv~=0);
    
    % Update the tract indexing to correspond with toi dimensions.
    hv = hv(find(hv~=0));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(toi(find(group == 1), :));
    m_wtwg2 = nanmean(toi(find(group == 2), :));
    m_wtwg3 = nanmean(toi(find(group == 3), :));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(toi(find(group == 1), :));
    std_wtwg2 = nanstd(toi(find(group == 2), :));
    std_wtwg3 = nanstd(toi(find(group == 3), :));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2;
    r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3;
    r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, [size(toi(find(group==1)))]), repmat(r_max_wtwg2, size(toi(find(group==2)))), ...
        repmat(r_max_wtwg3, size(toi(find(group==3)))));
    r_min = cat(1, repmat(r_min_wtwg1, [size(toi(find(group==1)))]), repmat(r_min_wtwg2, size(toi(find(group==2)))), ...
        repmat(r_min_wtwg3, size(toi(find(group==3)))));
    
    % Replace outliers with NaN.
    disp([wm_measure_here{w}]);
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
    t_out = array2table(cat(2, sub, group, cov_age, toi(:, 1:end), nanmean(toi(:, hv == 1), 2), nanmean(toi(:, hv == 2), 2)), 'VariableNames', ...
        {'subID', 'Age_group', 'Age_mo', list_tract{:}, 'meanH', 'meanV'});
    
    % Write.
    writetable(t_out, [rootDir 'LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '_test32104.csv']);
    
    % Output z-scored file for SPSS.
    toi_z = (nanmean(toi, 1) - toi)./nanstd(toi, [], 1);
    temp = nanmean(toi(:, hv == 1), 2);
    toi_meanh_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    temp = nanmean(toi(:, hv == 2), 2);
    toi_meanv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    t_out_z = array2table(cat(2, sub, group, cov_age, toi_z(:, 1:end), toi_meanh_z, toi_meanv_z), 'VariableNames', ...
        {'subID', 'Age_group', 'Age_mo', list_tract{:}, 'meanH', 'meanV'});
    writetable(t_out_z, [rootDir 'LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '_test32104_z.csv']);
    
end


