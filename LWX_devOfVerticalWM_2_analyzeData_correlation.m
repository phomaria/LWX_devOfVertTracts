% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

wm_measure = 'md'; %fa, od, icvf, isovf
y_min = .7; y_max = 1.15;

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
% addpath(genpath([rootDir 'proj-5a74d1e26ed91402ce400cca/']));

beh_measure = 'age'; %age, lit, vm, fm

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_tractz.mat'])
p_val = 0.05; p_bonf = p_val/length(list_tract);

save_figures = 'yes';

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

% Select the measurements of the tracts that I care about and categorize them into h or v. Convert all zeros to NaN.
ht = wm_childrenOnly(:, h_idx); ht(ht==0) = NaN;
vt = wm_childrenOnly(:, v_idx); vt(vt==0) = NaN;

% group mean z-score: performed within-category, across-subjects to control for subject-level
% differences in 'reactivity' while keeping subject-level differences in the WM measurement of interest
% ht(isnan(ht))=0; % set NaNs to the mean for now
z_ht = (nanmean(ht, 1) - ht)./nanstd(ht, [], 1);

% vt(isnan(vt))=0; % set NaNs to the mean for now
z_vt = (nanmean(vt, 1) - vt)./nanstd(vt, [], 1);

% Subset list_tracts so that we can call the correct tract names later.
list_tract_ht = list_tract(h_idx, :);
list_tract_vt = list_tract(v_idx, :);

c = colorcube;
%% HORIZONTAL TRACTS
% Does FA correlate with age in children?
count = 0;
figure(1); r_out = zeros(size(list_tract_ht)); p_out = zeros(size(list_tract_ht));
% Go through each tract, one at a time.
for t = 1:length(list_tract_ht)
    
    count = count + 4; 
        
    % Plot only data points that would have been included in the corelation. 
    scatter(cov_age_childrenOnly, ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
    if t == 1; hold on; end
    xlim([min(cov_age_childrenOnly) min(cov_age_childrenOnly)+10])
    if strcmp(beh_measure, 'age')
        [r, p] = plotcorr2(cov_age_childrenOnly, ht(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure, ', zscored'], [], list_tract_ht{t}, c(count, :));
    else
        [r, p] = plotcorr2(cov_age_childrenOnly, ht(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_ht{t}, c(count, :));
    end
    r_out(t) = r;
    p_out(t) = p;

end

title('Horizontal Tracts')
ylim([y_min y_max]);

if strcmp(save_figures, 'yes')
    
    print([rootDir 'plots/plot_correlation_age_' wm_measure '_horizontaltracts'], '-dpng')
    print([rootDir 'plots/eps/plot_correlation_age_' wm_measure '_horizontaltracts'], '-depsc')
    
end
hold off
disp(list_tract_ht)
disp(r_out)
disp(p_out)

clear r_out p_out

%% VERTICAL TRACTS
% Does FA correlate with age in children?
count = 64;
figure(2); r_out = zeros(size(list_tract_vt)); p_out = zeros(size(list_tract_vt));
% Go through each tract, one at a time.
for t = 1:length(list_tract_vt)
    
    count = count - 3; 
        
    % Plot only data points that would have been included in the corelation. 
    scatter(cov_age_childrenOnly, vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
    if t == 1; hold on; end
    xlim([min(cov_age_childrenOnly) min(cov_age_childrenOnly)+10])
    if strcmp(beh_measure, 'age')
        [r, p] = plotcorr2(cov_age_childrenOnly, vt(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure, ', zscored'], [], list_tract_vt{t}, c(count, :));
    else
        [r, p] = plotcorr2(cov_age_childrenOnly, vt(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_vt{t}, c(count, :));
    end
    r_out(t) = r;
    p_out(t) = p;
end

title('Vertical Tracts')
ylim([y_min y_max]);

disp(list_tract_ht)
disp(r_out)
disp(p_out)

if strcmp(save_figures, 'yes')
    
    print([rootDir 'plots/plot_correlation_age_' wm_measure '_verticaltracts'], '-dpng')
    print([rootDir 'plots/eps/plot_correlation_age_' wm_measure '_verticaltracts'], '-depsc')
    
end
hold off