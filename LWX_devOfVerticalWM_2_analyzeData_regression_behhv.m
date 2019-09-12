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
% addpath(genpath([rootDir 'proj-5a74d1e26ed91402ce400cca/']));

wm_measure = 'fa'; %fa, ad, rd, md, od, icvf, isovf
beh_measure = {'age', 'lit', 'vm', 'fm'}; %NOTE: this only uses age right now.

%% Tractography

% Read in data.
    load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure '_raw.mat']))
    
    % Convert into array and header for ease.
    data_all_in = table2array(data_tbl);
    data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Get grouping variable. NOTE: need to add lit, vm, and fm.
%     if strcmp(beh_measure{1}, 'age')
        
        group = data_tbl.gp_age;
        
%     end

% Get index matrices for hypothesis-driven grouping of WM tracts. 
for k = 1:length(data_all_in_header)
        
        % Indices of horizontal tracts.
        h_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
            || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF') ...
            || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
            || strcmp(data_all_in_header{k}, 'leftArc') || strcmp(data_all_in_header{k}, 'rightArc') ...
            || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3');
        
        % Indices of vertical tracts.
        v_idx(k) = strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
            || strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
            || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
            || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
            || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang');
        
end

% SELECT the measurements of the tracts that I care (h and v) and 
% the subjects that I care about (non-adults). 
% Categorize into h or v. Convert all zeros to NaN.
ht = data_all_in(group ~= 3, h_idx); ht(ht==0) = NaN;
vt = data_all_in(group ~= 3, v_idx); vt(vt==0) = NaN;

% AVERAGE the WM propery of interest for each subject averaged across tracts within horizontal or vertical categories.
categorymean_h = nanmean(ht, 2);
categorymean_v = nanmean(vt, 2);

% ZSCORE the averages based on tract mean and standard deviation (i.e., dim=1) using the sample standard deviation (i.e., flag=0)
z_categorymean_h = zscore(categorymean_h, 0, 1);
z_categorymean_v = zscore(categorymean_v, 0, 1);

% Order the WM predictor to conform to the required data structure for LMM: https://www.mathworks.com/help/stats/prepare-data-for-linear-mixed-effects-models.html
temp = transpose(cat(2, z_categorymean_h, z_categorymean_v));
z_categorymean = temp(:); clear temp

%% Behavior.

% Subselect only data for children and for the covariates of interest. 
beh = cat(2, data_tbl.cov_age(group~=3), data_tbl.c_lit(group~=3), data_tbl.c_vm(group~=3), data_tbl.c_fm(group~=3), ...
    data_tbl.cov_sex(group~=3));

% Get measure-specific z-scores (even for age).
z_beh = (beh - nanmean(beh))./nanstd(beh);

%% Perform Multiple Linear Regression: [WMh WMv] ~ age + lit + vm + fm + 1|Subject.

% RESPONSE VARIABLE: Define behavior the response variable.
y = cat(1, z_categorymean_h, z_categorymean_v);

% TEST: Add noise to respons variable.
% y = y + rand(size(y));

% FIXED EFFECTS: define design matrix.
X = blkdiag(z_beh, z_beh);

% RANDOM EFFECT: define random intercept.
Z = ones(size(y));

% GROUPING VARIABLE: subject.
G = repmat(data_tbl.subID(group ~= 3), [2 1]); 

% Convert to table for LME.
tbl = array2table(cat(2, y, X, Z, G), 'VariableNames', {'wm', 'age_h', 'lit_h', 'vm_h', 'fm_h', 'sex_h', 'age_v', 'lit_v', 'vm_v', 'fm_v', 'sex_v', 'constant', 'sub'});

% Perform fitting procedure with fitlme: 'Recog ~ WMh + Wmv + 1|subj';
% fitlemematrix uses fitlme but does not require table format
mdl = fitlmematrix([ones(size(y)) X], y, Z, G, 'FixedEffectPredictors', {'intercept', 'age_Horizontal', 'lit_Horizontal', ...
    'vm_Horizontal', 'fm_Horizontal', 'sex_Horizontal', 'age_Vertical', 'lit_Vertical',  'vm_Vertical',  'fm_Vertical', 'sex_Vertical'}, ...
    'RandomEffectPredictors', {'intercept'}, 'RandomEffectGroups', {'sub'}) 

tag = 1;

format short

% Display R2.
disp(['Rsquared is: ' num2str(mdl.Rsquared.Ordinary) '.']);

% Display f- and p-statistics for the overall model.
[p3,f3] = coefTest(mdl);
disp(['Model: F = ' num2str(f3) ', p = ' num2str(p3)])

% Display betas.
disp(mdl.Coefficients)

% Posthoc comparison for v > h.
[pVal, F, df1, df2] = coefTest(mdl, [0 -1 -1 -1 -1 0 1 1 1 1 0]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for All Behavioral Measures: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Posthoc comparison for v > h, age.
[pVal, F, df1, df2] = coefTest(mdl, [0 -1 0 0 0 0 1 0 0 0 0]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for Age: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Posthoc comparison for v > h, literacy.
[pVal, F, df1, df2] = coefTest(mdl, [0 0 -1 0 0 0 0 1 0 0 0]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for Literacy: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Posthoc comparison for v > h, visual-motor.
[pVal, F, df1, df2] = coefTest(mdl, [0 0 0 -1 0 0 0 0 1 0 0]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for Visual-Motor Skill: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Posthoc comparison for v > h, fine-motor.
[pVal, F, df1, df2] = coefTest(mdl, [0 0 0 0 -1 0 0 0 0 1 0]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for Fine-Motor SKill: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Posthoc comparison for v > h, sex.
[pVal, F, df1, df2] = coefTest(mdl, [0 0 0 0 0 -1 0 0 0 0 1]);
disp(['Planned comparison All Vertical Tracts > All Horizontal Tracts for Sex: F(' ...
    num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', p = ' num2str(pVal) '.'])

% Set significance thresholds for post-hoc comparisons.
p_val = 0.05; p_bonf = p_val/5;
disp(['Bonferroni-corrected p-value: ' num2str(p_bonf)])

%% Visualize.

figure(1)
xticks = linspace(0, 1, length(mdl.Coefficients.Estimate));
markersize = 8; buffer = 0.2; fontsize = 14; linesize = 3;
% y_min = min(mdl.Coefficients.SE) - 1; y_max = max(mdl.Coefficients.SE) + 1;
y_min = -1.5; y_max = 1.5;

% Age
h1 = plot(xticks([1, 6]), mdl.Coefficients.Estimate([tag+1, tag+6]), 'ko', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize-1, 'HandleVisibility', 'off'); hold on;
errorbars = errorbar(xticks([1, 6]), mdl.Coefficients.Estimate([tag+1, tag+6]), mdl.Coefficients.SE([tag+1, tag+6]));
set(errorbars, 'LineStyle', 'none');
set(errorbars,  'LineWidth', 1, 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);

% Literacy
h2 = plot(xticks([2, 7]), mdl.Coefficients.Estimate([tag+2, tag+7]), 'square', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);
errorbars = errorbar(xticks([2, 7]), mdl.Coefficients.Estimate([tag+2, tag+7]), mdl.Coefficients.SE([tag+2, tag+7]));
set(errorbars, 'LineStyle', 'none');
set(errorbars,  'LineWidth', 1, 'Color', 'k', 'Marker', 'square', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);

% Visual-motor
h3 = plot(xticks([3, 8]), mdl.Coefficients.Estimate([tag+3, tag+8]), '^', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);
errorbars = errorbar(xticks([3, 8]), mdl.Coefficients.Estimate([tag+3, tag+8]), mdl.Coefficients.SE([tag+3, tag+8]));
set(errorbars, 'LineStyle', 'none');
set(errorbars,  'LineWidth', 1, 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);

% Fine-motor
h4 = plot(xticks([4, 9]), mdl.Coefficients.Estimate([tag+4, tag+9]), 'diamond', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);
errorbars = errorbar(xticks([4, 9]), mdl.Coefficients.Estimate([tag+4, tag+9]), mdl.Coefficients.SE([tag+9, tag+9]));
set(errorbars, 'LineStyle', 'none');
set(errorbars,  'LineWidth', 1, 'Color', 'k', 'Marker', 'diamond', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', markersize);

% Dividing lines.
if tag == 1 % random intercept for subject was not used
    
    plot([(max(xticks)+3*buffer)*0.27 (max(xticks)+3*buffer)*0.27], [y_min y_max], ':k', 'LineWidth', linesize)
    xlim([min(xticks)-buffer max(xticks)]); ylim([y_min y_max]);
    plot(linspace(min(xticks)-buffer, max(xticks), 10), zeros(1, 10), 'k-')
    text(0, y_min+buffer, 'Horizontal', 'FontSize', fontsize); text(0.60, y_min+buffer, 'Vertical', 'FontSize', fontsize);
    
elseif tag == 0 % random intercept for subject was used
    
    plot([(max(xticks)+4*buffer)*0.29 (max(xticks)+4*buffer)*0.29], [y_min-buffer y_max], ':k', 'LineWidth', linesize)
    xlim([min(xticks)-buffer max(xticks)+buffer]); ylim([y_min-2*buffer y_max+0.5*buffer]);
    plot(linspace(min(xticks)-buffer, max(xticks)+buffer, 10), zeros(1, 10), 'k-', 'LineWidth', linesize)
    text(0.05, y_min-buffer, 'Horizontal', 'FontSize', fontsize); text(0.80, y_min-buffer, 'Vertical', 'FontSize', fontsize);
    
end

legend([h1 h2 h3 h4], beh_measure, 'Location', 'southoutside', 'Orientation', 'horizontal');
set(gca,'xtick',[]); set(gca,'xticklabel',{[]}); set(gca,'FontSize', fontsize)
ylabel('Beta Estimates +/- SE');
box off;

if strcmp(wm_measure, 'fa')
    fig_title = 'Fractional Anisotropy';
elseif strcmp(wm_measure, 'od')
    fig_title = 'Orientation Dispersion';
elseif strcmp(wm_measure, 'icvf')
    fig_title = 'Neurite Density (ICVF)';
elseif strcmp(wm_measure, 'isovf')
    fig_title = 'Isotropic Volume Fraction (ISOVF)';
elseif strcmp(wm_measure, 'ad')
    fig_title = 'Axial Diffusivity (AD)';
elseif strcmp(wm_measure, 'md')
    fig_title = 'Mean Diffusivity (MD)';
elseif strcmp(wm_measure, 'rd')
    fig_title = 'Radial Diffusivity (RD)';
end
title([fig_title ', r2 = ' num2str(mdl.Rsquared.Ordinary)]);

print(fullfile(rootDir, 'plots', ['plot_mdlEstimates_' wm_measure '_hv']), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_mdlEstimates_' wm_measure '_hv']), '-depsc')

hold off;

