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

wm_measure = 'fa'; %fa [.4, .7], ad [1, 2], rd [.4, .8], md [.6, 1.1], icvf [.4, 1], isovf [0, .3], od [.1, .4]
y_min = 0.4; y_max = 0.7;
save_figures = 'yes';

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
% addpath(genpath([rootDir 'proj-5a74d1e26ed91402ce400cca/']));

beh_measure = 'age'; %age, lit, vm, fm

% Read in data.
load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure '_raw.mat']));

% Convert into array and header for ease.
data_all_in = table2array(data_tbl);
data_all_in_header = data_tbl.Properties.VariableNames;

% Get grouping variable. NOTE: need to add lit, vm, and fm. This means also
% that measure_childrenOnly_z needs to be defined for the lit, vm, and fm cases.
if strcmp(beh_measure, 'age')
    
    group = data_tbl.gp_age;
    
end

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

% Select the measurements of the tracts that I care (h and v) and 
% the subjects that I care about (non-adults). 
% Categorize into h or v. Convert all zeros to NaN.
ht = data_all_in(group ~= 3, h_idx); ht(ht==0) = NaN;
vt = data_all_in(group ~= 3, v_idx); vt(vt==0) = NaN;

% Do the same for the adults just to show it on the plot; only makes sense for age.
if strcmp(beh_measure, 'age')
    
    ht_adult = data_all_in(group == 3, h_idx); ht_adult(ht_adult==0) = NaN;
    vt_adult = data_all_in(group == 3, v_idx); vt_adult(vt_adult==0) = NaN;
    
end

% group mean z-score: performed within-category, across-subjects to control for subject-level
% differences in 'reactivity' while keeping subject-level differences in the WM measurement of interest
% ht(isnan(ht))=0; % set NaNs to the mean for now
z_ht = (nanmean(ht, 1) - ht)./nanstd(ht, [], 1);

% vt(isnan(vt))=0; % set NaNs to the mean for now
z_vt = (nanmean(vt, 1) - vt)./nanstd(vt, [], 1);

% Subset list_tracts so that we can call the correct tract names later.
list_tract_ht = data_all_in_header(h_idx);
list_tract_vt = data_all_in_header(v_idx);

% Get z-scored age.
z_cov_age_childrenOnly = (nanmean(data_tbl.cov_age(group~=3), 1) - data_tbl.cov_age(group~=3))./nanstd(data_tbl.cov_age(group~=3), [], 1);

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
        scatter(data_tbl.cov_age(group~=3), ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(ht_adult(:, t))), ht_adult(:, t), 'MarkerEdgeColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(data_tbl.cov_age(group~=3), ht(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_ht{t}, c(count, :), 5);
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
        scatter(data_tbl.cov_age(group~=3), vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(vt_adult(:, t))), vt_adult(:, t), 'MarkerEdgeColor', c(count, :))
        if t == 1; hold on; end
        [r, p, ~, slope] = plotcorr2(data_tbl.cov_age(group~=3), vt(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_vt{t}, c(count, :), 5);
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
        scatter(data_tbl.cov_age(group~=3), ht(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(ht_adult(:, t))), ht_adult(:, t), 'MarkerEdgeColor', c(count, :))
        [r, p, statstr] = plotcorr2(data_tbl.cov_age(group~=3), ht(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_ht{t}, c(count, :), 5);
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
        scatter(data_tbl.cov_age(group~=3), vt(:, t), 'filled', 'MarkerFaceColor', c(count, :))
        scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(vt_adult(:, t))), vt_adult(:, t), 'MarkerEdgeColor', c(count, :))
        [r, p, statstr] = plotcorr2(data_tbl.cov_age(group~=3), vt(:, t), [beh_measure ' (mo)'], ['Average ' wm_measure], [], list_tract_vt{t}, c(count, :), 5);
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
    scatter(repmat(data_tbl.cov_age(group~=3), [size(ht, 2) 1]), ht(:), 'filled', 'MarkerFaceColor', c(count, :))
    hold on;
    scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(ht_adult(:))), ht_adult(:), 'MarkerEdgeColor', c(count, :))
    [r, p, ~, slope] = plotcorr2(repmat(data_tbl.cov_age(group~=3), [size(ht, 2) 1]), ht(:), [beh_measure ' (mo)'], ['Average ' wm_measure], [], 'Horizontal', c(count, :), 5);
    ylim([y_min y_max]);
    title(['Horizontal Tracts, slope = ' num2str(slope) ', r = ' num2str(r) ', p = ' num2str(p)]);
    hold off;
    
    % Vertical
    figure(x+1)
    scatter(repmat(data_tbl.cov_age(group~=3), [size(vt, 2) 1]), vt(:), 'filled', 'MarkerFaceColor', c(count+30, :))
    hold on;
    scatter(repmat(max(data_tbl.cov_age(group~=3) + 4), size(vt_adult(:))), vt_adult(:), 'MarkerEdgeColor', c(count+30, :))
    [r, p, ~, slope] = plotcorr2(repmat(data_tbl.cov_age(group~=3), [size(vt, 2) 1]), vt(:), [beh_measure ' (mo)'], ['Average ' wm_measure], [], 'Vertical', c(count+30, :), 5);
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
