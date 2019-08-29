% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

control_age = 'yes';

wm_measure = 'fa'; %fa, od, icvf, isovf
y_min = .4; y_max = .7;

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

if strcmp(beh_measure, 'age')
    
    ht_adult = wm(find(group == 3), h_idx); ht_adult(ht_adult==0) = NaN;
    vt_adult = wm(find(group == 3), v_idx); vt_adult(vt_adult==0) = NaN;
    
end

% group mean z-score: performed within-category, across-subjects to control for subject-level
% differences in 'reactivity' while keeping subject-level differences in the WM measurement of interest
% ht(isnan(ht))=0; % set NaNs to the mean for now
z_ht = (nanmean(ht, 1) - ht)./nanstd(ht, [], 1);

% vt(isnan(vt))=0; % set NaNs to the mean for now
z_vt = (nanmean(vt, 1) - vt)./nanstd(vt, [], 1);

% Subset list_tracts so that we can call the correct tract names later.
list_tract_ht = list_tract(h_idx, :);
list_tract_vt = list_tract(v_idx, :);

% Get z-scored age.
z_cov_age_childrenOnly = (nanmean(cov_age_childrenOnly, 1) - cov_age_childrenOnly)./nanstd(cov_age_childrenOnly, [], 1);

c = colorcube;
%% HORIZONTAL TRACTS
% Does FA correlate with age in children?
count = 0;
figure(1); r_out = zeros(size(list_tract_ht)); p_out = zeros(size(list_tract_ht));
% Go through each tract, one at a time.
for t = 1:length(list_tract_ht)
    
    count = count + 4;
    
    if strcmp(beh_measure, 'age')
        
        % Plot only data points that would have been included in the correlation.
        scatter(cov_age_childrenOnly, ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(cov_age_childrenOnly + 4), size(ht_adult(:, t))), ht_adult(:, t), 'MarkerEdgeColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(cov_age_childrenOnly, ht(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_ht{t}, c(count, :), 5);
        ylim([y_min y_max]);
        
    elseif strcmp(control_age, 'yes')
        
        % Plot only data points that would have been included in the correlation.
        scatter(measure_childrenOnly_z, z_ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotpartialcorr(measure_childrenOnly_z, z_ht(:, t), z_cov_age_childrenOnly, beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_ht{t}, c(count, :));
        xlim([-2 2]);
        ylim([-2 2]);
        
    else
        
        % Plot only data points that would have been included in the correlation.
        scatter(measure_childrenOnly_z, z_ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(measure_childrenOnly_z, z_ht(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_ht{t}, c(count, :), .2);
        xlim([-2 2]);
        ylim([-2 2]);
        
    end
    
    r_out(t) = r;
    p_out(t) = p;
    slope_out(t) = slope;
    
end

title('Horizontal Tracts')

if strcmp(save_figures, 'yes')
    
    print([rootDir 'plots/plot_correlation_' beh_measure '_' wm_measure '_horizontaltracts'], '-dpng')
    print([rootDir 'plots/eps/plot_correlation_' beh_measure '_' wm_measure '_horizontaltracts'], '-depsc')
    
end
hold off
disp(list_tract_ht)
disp(r_out)
disp(p_out)
slope_h = slope_out;

clear r_out p_out slope_out

%% VERTICAL TRACTS
% Does FA correlate with age in children?
count = 64;
figure(2); r_out = zeros(size(list_tract_vt)); p_out = zeros(size(list_tract_vt));
% Go through each tract, one at a time.
for t = 1:length(list_tract_vt)
    
    count = count - 3;
    
    if strcmp(beh_measure, 'age')
        
        % Plot only data points that would have been included in the corelation.
        scatter(cov_age_childrenOnly, vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(cov_age_childrenOnly + 4), size(vt_adult(:, t))), vt_adult(:, t), 'MarkerEdgeColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(cov_age_childrenOnly, vt(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_vt{t}, c(count, :), 5);
        ylim([y_min y_max]);
        
    elseif strcmp(control_age, 'yes')
        
        
        % Plot only data points that would have been included in the correlation.
        scatter(measure_childrenOnly_z, z_vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        if t == 1; hold on; end
        [r, p] = plotpartialcorr(measure_childrenOnly_z, z_vt(:, t), z_cov_age_childrenOnly, beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_vt{t}, c(count, :));
        xlim([-2 2]);
        ylim([-2 2]);
        
    else
        
        % Plot only data points that would have been included in the corelation.
        scatter(measure_childrenOnly_z, z_vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(measure_childrenOnly_z, z_vt(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_vt{t}, c(count, :), .2);
        xlim([-2 2]);
        ylim([-2 2]);
        
    end
    
    r_out(t) = r;
    p_out(t) = p;
    slope_out(t) = slope;
    
end

title('Vertical Tracts')

disp(list_tract_vt)
disp(r_out)
disp(p_out)
slope_v = slope_out;

if strcmp(save_figures, 'yes')
    
    print([rootDir 'plots/plot_correlation_' beh_measure '_' wm_measure '_verticaltracts'], '-dpng')
    print([rootDir 'plots/eps/plot_correlation_' beh_measure '_' wm_measure '_verticaltracts'], '-depsc')
    
end
hold off

%% EACH TRACT ONE AT A TIME
% Does FA correlate with age in children?
% Go through each tract, one at a time.
count = 0;
for t = 1:length(list_tract_ht)
    
    count = count + 4;
    
    figure(t+2); hold on;
   
    if strcmp(beh_measure, 'age')
        
        % Plot only data points that would have been included in the corelation.
        scatter(cov_age_childrenOnly, ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(cov_age_childrenOnly + 4), size(ht_adult(:, t))), ht_adult(:, t), 'MarkerEdgeColor', c(count, :))
        [r, p, statstr] = plotcorr2(cov_age_childrenOnly, ht(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_ht{t}, c(count, :), 5);
        ylim([y_min y_max]);
        
    elseif strcmp(control_age, 'yes')
        
        % Plot only data points that would have been included in the corelation.
        scatter(measure_childrenOnly_z, z_ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        [r, p, statstr] = plotpartialcorr(measure_childrenOnly_z, z_ht(:, t), z_cov_age_childrenOnly, beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_ht{t}, c(count, :));
        xlim([-2 2]);
        ylim([-2 2]);
        
    else
        
        % Plot only data points that would have been included in the corelation.
        scatter(measure_childrenOnly_z, z_ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        [r, p, statstr] = plotcorr2(measure_childrenOnly_z, z_ht(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_ht{t}, c(count, :), .2);
        xlim([-2 2]);
        ylim([-2 2]);
        
    end
    
    title([list_tract_ht{t} ', ' statstr])
    hold off
    
end

count = 64; x=t+2;
for t = 1:length(list_tract_vt)
    
    count = count - 3;
    
    figure(t+x); hold on;
    
    if strcmp(beh_measure, 'age')
        
        % Plot only data points that would have been included in the corelation.
        scatter(cov_age_childrenOnly, vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(cov_age_childrenOnly + 4), size(vt_adult(:, t))), vt_adult(:, t), 'MarkerEdgeColor', c(count, :))
        [r, p, statstr] = plotcorr2(cov_age_childrenOnly, vt(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_vt{t}, c(count, :), 5);
        ylim([y_min y_max]);
        
    elseif strcmp(control_age, 'yes')
        
        % Plot only data points that would have been included in the corelation.
        scatter(measure_childrenOnly_z, z_vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        [r, p, statstr] = plotpartialcorr(measure_childrenOnly_z, z_vt(:, t), z_cov_age_childrenOnly, beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_vt{t}, c(count, :));
        xlim([-2 2]);
        ylim([-2 2]);
        
    else
        
        % Plot only data points that would have been included in the corelation.
        scatter(measure_childrenOnly_z, z_vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        [r, p, statstr] = plotcorr2(measure_childrenOnly_z, z_vt(:, t), beh_measure, ['Average ' wm_measure, ', z-scored'], [], list_tract_vt{t}, c(count, :), .2);
        xlim([-2 2]);
        ylim([-2 2]);
        
    end
    
    title([list_tract_vt{t} ', ' statstr])
    hold off;
    
end

x = t+2;
% Plot just vertical correlation and just horizontal correlation on same plot.
if strcmp(beh_measure, 'age')
    
    % Horizontal
    figure(x)
    scatter(repmat(cov_age_childrenOnly, [size(ht, 2) 1]), ht(:), 'filled', 'MarkerFaceColor', c(count, :))
    hold on;
    scatter(repmat(max(cov_age_childrenOnly + 4), size(ht_adult(:))), ht_adult(:), 'MarkerEdgeColor', c(count, :))
    [r, p, ~, slope] = plotcorr2(repmat(cov_age_childrenOnly, [size(ht, 2) 1]), ht(:), [beh_measure ' (mo)'], ['Average ' wm_measure], [], 'Horizontal', c(count, :), 5);
    title(['Horizontal Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
    % Vertical
    figure(x+1)
    scatter(repmat(cov_age_childrenOnly, [size(vt, 2) 1]), vt(:), 'filled', 'MarkerFaceColor', c(count+30, :))
    hold on;
    scatter(repmat(max(cov_age_childrenOnly + 4), size(vt_adult(:))), vt_adult(:), 'MarkerEdgeColor', c(count+30, :))
    [r, p, ~, slope] = plotcorr2(repmat(cov_age_childrenOnly, [size(vt, 2) 1]), vt(:), [beh_measure ' (mo)'], ['Average ' wm_measure], [], 'Vertical', c(count+30, :), 5);
    ylim([y_min y_max]);
    title(['Vertical Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
elseif strcmp(control_age, 'yes')

    % Horizontal
    figure(x)
    scatter(repmat(measure_childrenOnly_z, [size(z_ht, 2) 1]), z_ht(:), 'filled', 'MarkerFaceColor', c(count, :))
    hold on;
    [r, p, ~, slope] = plotpartialcorr(repmat(measure_childrenOnly_z, [size(z_ht, 2) 1]), z_ht(:), repmat(z_cov_age_childrenOnly, [size(ht, 2) 1]), beh_measure, ['Average ' wm_measure, ', z-scored'], [], 'Horizontal', c(count, :));
    xlim([-2 2]);
    ylim([-2 2]);
    title(['Horizontal Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
    % Vertical
    figure(x+1)
    scatter(repmat(measure_childrenOnly_z, [size(z_vt, 2) 1]), z_vt(:), 'filled', 'MarkerFaceColor', c(count+30, :))
    hold on;
    [r, p, ~, slope] = plotpartialcorr(repmat(measure_childrenOnly_z, [size(z_vt, 2) 1]), z_vt(:), repmat(z_cov_age_childrenOnly, [size(vt, 2) 1]), beh_measure, ['Average ' wm_measure, ', z-scored'], [], 'Vertical', c(count+30, :));
    xlim([-2 2]);
    ylim([-2 2]);
    title(['Vertical Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
else
    
    % Horizontal
    figure(x)
    scatter(measure_childrenOnly_z, z_ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
    hold on;
    [r, p, ~, slope] = plotcorr2(repmat(measure_childrenOnly_z, [size(z_ht, 2) 1]), z_ht(:), beh_measure, ['Average ' wm_measure, ', z-scored'], [], 'Horizontal', c(count, :));
    xlim([-2 2]);
    ylim([-2 2]);
    title(['Horizontal Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
    % Vertical
    figure(x+1)
    scatter(measure_childrenOnly_z, z_vt(:, t), 'filled', 'MarkerFaceColor', c(count+30, :))
    hold on;
    [r, p, ~, slope] = plotcorr2(repmat(measure_childrenOnly_z, [size(z_vt, 2) 1]), z_vt(:), beh_measure, ['Average ' wm_measure, ', z-scored'], [], 'Vertical', c(count+30, :));
    xlim([-2 2]);
    ylim([-2 2]);
    title(['Vertical Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
end

% NOTE: need to change this to paired sample ttest with unequal samples.
[h,p,ci,stats] = ttest2(slope_h, slope_v)
