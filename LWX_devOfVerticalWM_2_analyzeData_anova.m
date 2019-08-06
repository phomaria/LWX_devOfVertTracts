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
save_figures = 'no';

beh_measure = 'age'; %age, lit, vm, fm
wm_measure = 'fa'; %fa, ad, md, rd, od, icvf, isovf

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_tractz_test32104.mat'])
clearvars -except save_figures wm group list_tract wm_measure rootDir covariates sub beh_measure beh_measures

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

% Reorganize data to prepare for two-way repeated measures anova using anovan: dependent variable.
toi_big = toi(:);

% Reorganize data to prepare for two-way repeated measures anova using anovan: grouping variables (independent variables).
gp_hv = [];
for r = 1:size(toi, 2)
    
    gp_temp = repmat(hv(r), [size(toi, 1), 1]);
    
    gp_hv = cat(1, gp_hv, gp_temp);
    
end
clear gp_temp
gp_sub = repmat(sub, [size(toi, 2), 1]);
gp_age = repmat(group, [size(toi, 2), 1]);

%% ==================== FIRST ANOVA ==================== %%

% --- 3(gp_age = 1 (youngchild), 2 (old child), 3 (adult)) X 2(gp_hv = 1 (horizontal), 2 (vertical)) Repeated Measures ANOVA using anovan
% Note: anovan for RM requires the subID to be coded, noted as a random
% variable and then for the subID to be nested under the subject grouping
% variable (gp_age).
% Test the 2-ways interactions.

[p, tbl, stats] = anovan(toi_big, {gp_age, gp_hv, gp_sub}, 'random', 3, 'nested', [0 0 0 ; 0 0 0 ; 1 0 0], ...
    'model','interaction', 'varnames', {'Age Group', 'Tract Orientation', 'subID'}); %, 'display', 'off');

% [p, tbl, stats] = anovan(toi_big, {gp_age, gp_hv}, ...
%     'model','interaction', 'varnames', {'Age Group', 'Tract Orientation'}); %, 'display', 'off');

% toi_h = nanmean(toi(:, find(hv==1)), 2); 
% toi_v = nanmean(toi(:, find(hv==2)), 2); 
% bt_tbl = table(cat(1, toi_h, toi_v), gp_age(1:size(toi, 1)*2), gp_hv(1:size(toi, 1)*2), 'VariableNames', {'WM', 'AgeGroup', 'TractOrientation'});
% rm = fitrm(bt_tbl, 'WM~AgeGroup*TractOrientation', 'WithinDesign', [1 2])
% ranovatbl = ranova(rm, 'WithinModel', 'w1+w2')
% anova(rm)

% toi_hv(:, 1) = nanmean(toi(:, find(hv==1)), 2); 
% toi_hv(:, 2) = nanmean(toi(:, find(hv==2)), 2); 
% % tbl = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
% tbl = simple_mixed_anova(toi_hv, group, {'AgeGroup'}, {'TractOrientation'})



% % Looking for an interaction between Age Group and Tract Orientation,
% % so report on p(4), in this case.
% p_indx = 4;
% disp(p)

% Visualize.
figure(1)
c = colorcube;
yc_color = c(1:11, :); 
oc_color = c(12:23, :); 
a_color = c(24:35, :);
G = findgroups(gp_sub); 
gradient = transpose(linspace(.5, 1, length(unique(G(find(gp_age==1 & gp_hv==1))))));
% Visualize: young children, horizontal
h1 = gscatter(ones(size(find(gp_age==1 & gp_hv==1), 1), 1), toi_big(find(gp_age==1 & gp_hv==1)), G(find(gp_age==1 & gp_hv==1)), yc_color);
%'.', 'MarkerFaceColor', yc_color, 'MarkerEdgeColor', yc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
hold on;
plot(linspace(.8, 1.2, 5), repmat(nanmean(toi_big(find(gp_age==1 & gp_hv==1))), [1 5]), 'k-')

% Visualize: older children, horizontal
h2 = gscatter(2*ones(size(find(gp_age==2 & gp_hv==1), 1), 1), toi_big(find(gp_age==2 & gp_hv==1)), G(find(gp_age==2 & gp_hv==1)), oc_color);
%'.', 'MarkerFaceColor', oc_color, 'MarkerEdgeColor', oc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
plot(linspace(1.8, 2.2, 5), repmat(nanmean(toi_big(find(gp_age==2 & gp_hv==1))), [1 5]), 'k-')

% Visualize: adults, horizontal
h3 = gscatter(3*ones(size(find(gp_age==3 & gp_hv==1), 1), 1), toi_big(find(gp_age==3 & gp_hv==1)), G(find(gp_age==3 & gp_hv==1)), a_color);
%'.', 'MarkerFaceColor', a_color, 'MarkerEdgeColor', a_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
plot(linspace(2.8, 3.2, 5), repmat(nanmean(toi_big(find(gp_age==3 & gp_hv==1))), [1 5]), 'k-')

 % Visualize: young children, vertical
gscatter(4*ones(size(find(gp_age==1 & gp_hv==2), 1), 1), toi_big(find(gp_age==1 & gp_hv==2)), G(find(gp_age==1 & gp_hv==2)), yc_color);
%'.', 'MarkerFaceColor', yc_color, 'MarkerEdgeColor', yc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
plot(linspace(3.8, 4.2, 5), repmat(nanmean(toi_big(find(gp_age==1 & gp_hv==2))), [1 5]), 'k-')

% Visualize: older children, vertical
gscatter(5*ones(size(find(gp_age==2 & gp_hv==2), 1), 1), toi_big(find(gp_age==2 & gp_hv==2)), G(find(gp_age==2 & gp_hv==2)), oc_color);
%'.', 'MarkerFaceColor', oc_color, 'MarkerEdgeColor', oc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
plot(linspace(4.8, 5.2, 5), repmat(nanmean(toi_big(find(gp_age==2 & gp_hv==2))), [1 5]), 'k-')

% Visualize: adults, vertical
gscatter(6*ones(size(find(gp_age==3 & gp_hv==2), 1), 1), toi_big(find(gp_age==3 & gp_hv==2)), G(find(gp_age==3 & gp_hv==2)), a_color);
%'.', 'MarkerFaceColor', a_color, 'MarkerEdgeColor', a_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
plot(linspace(5.8, 6.2, 5), repmat(nanmean(toi_big(find(gp_age==3 & gp_hv==2))), [1 5]), 'k-')
b_plot = gca; legend(b_plot,'off');

% Dividing lines.
buffer = 0.3; linesize = 2; fontsize= 10;
plot(repmat(3.5, [1 5]), linspace(0, 5, 5),  ':k', 'LineWidth', linesize)

set(gca,'xtick',[1 2 3 4 5 6]); 
% set(gca,'xticklabel',{'children', 'children', 'adults', 'children', 'children', 'adults'});...
%     '4.5-6.0 yrs.', '6.0-8.5 yrs.', ' ', '4.5-6.0 yrs.', '6.0-8.5 yrs.', ' '}); 
labels = {'Children 4.5-6.0yrs.','Children 6.0-8.5yrs.','Adults'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
a = gca;
a.XTickLabel = labels;
set(gca,'FontSize', fontsize);

if strcmp(wm_measure, 'fa')
    fig_title = 'Fractional Anisotropy';
    ylim([0.3 0.75]);
elseif strcmp(wm_measure, 'od')
    fig_title = 'Orientation Dispersion';
    ylim([0 0.5]);
elseif strcmp(wm_measure, 'icvf')
    fig_title = 'Neurite Density (ICVF)';
    ylim([0.5 1]);
elseif strcmp(wm_measure, 'isovf')
    fig_title = 'Isometric Volume Fraction (ISOVF)';
    ylim([0 .3]);
elseif strcmp(wm_measure, 'ad')
    fig_title = 'Axial Diffusivity (AD)';
    ylim([0.6 1.8]);
elseif strcmp(wm_measure, 'md')
    fig_title = 'Mean Diffusivity (MD)';
    ylim([0 1.3]);
elseif strcmp(wm_measure, 'rd')
    fig_title = 'Radial Diffusivity (RD)';
    ylim([0 0.9]);
end
xlim([0.5 6.5]);
text(1.5, min(yticks+.05), 'Horizontal', 'FontSize', fontsize); text(4.6, min(yticks+.05), 'Vertical', 'FontSize', fontsize);
title([fig_title])
box off;

print([rootDir 'plots/plot_3x2anova_' wm_measure '_hv'], '-dpng')
print([rootDir 'plots/eps/plot_3x2anova_' wm_measure '_hv'], '-depsc')

hold off;
mc = 10;

disp('=======================')
disp('Planned Comparisons')
disp(' ')

% Check for old children vs. adult.
% FOLLOW-UP: Is WM in h different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==2 & gp_hv==1)), toi_big(find(gp_age==3 & gp_hv==1)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) > nanmean(toi_big(find(gp_age==2 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is greater in adults than in older children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) < nanmean(toi_big(find(gp_age==2 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is lower in adults than in older children, p= ' num2str(p2) '.'])
else
    disp(['No difference in horizontal tracts between adults and older children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in v different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==2 & gp_hv==2)), toi_big(find(gp_age==3 & gp_hv==2)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==2))) > nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is greater in adults than in older children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==2))) < nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is lower in adults than in older children, p= ' num2str(p2) '.'])
else
    disp(['No difference in vertical tracts between adults and older children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in children different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age==2 & gp_hv==1)), toi_big(find(gp_age==2 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) > nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in older children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) < nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in older children, p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in older children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in adults different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age==3 & gp_hv==1)), toi_big(find(gp_age==3 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) > nanmean(toi_big(find(gp_age==3 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in adults, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) < nanmean(toi_big(find(gp_age==3 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in adults, p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in adults, p= ' num2str(p2) '.'])
end

disp(' ')

% FOLLOW-UP: Is WM in h different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==1 & gp_hv==1)), toi_big(find(gp_age==2 & gp_hv==1)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) > nanmean(toi_big(find(gp_age==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is greater in older than younger children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) < nanmean(toi_big(find(gp_age==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is lower in older than younger children, p= ' num2str(p2) '.'])
else
    disp(['No difference in horizontal tracts between older and younger children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in v different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==1 & gp_hv==2)), toi_big(find(gp_age==2 & gp_hv==2)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==2))) > nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is greater in older than younger children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==2))) < nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is lower in older than younger children, p= ' num2str(p2) '.'])
else
    disp(['No difference in vertical tracts between older and younger children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in children different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age==1 & gp_hv==1)), toi_big(find(gp_age==1 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==1 & gp_hv==1))) > nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in younger children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==1 & gp_hv==1))) < nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in younger children, p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in younger children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in adults different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age==2 & gp_hv==1)), toi_big(find(gp_age==2 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) > nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in older children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==2 & gp_hv==1))) < nanmean(toi_big(find(gp_age==2 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in older children, p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in older children, p= ' num2str(p2) '.'])
end

disp(' ')

% FOLLOW-UP: Is WM in h different in younger children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==3 & gp_hv==1)), toi_big(find(gp_age==1 & gp_hv==1)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) > nanmean(toi_big(find(gp_age==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is greater in adults than younger children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==1))) < nanmean(toi_big(find(gp_age==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is lower in adults than younger children, p= ' num2str(p2) '.'])
else
    disp(['No difference in horizontal tracts between adults and younger children, p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in v different in younger children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age==3 & gp_hv==2)), toi_big(find(gp_age==1 & gp_hv==2)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==2))) > nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is greater in adults than younger children, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age==3 & gp_hv==2))) < nanmean(toi_big(find(gp_age==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is lower in adults than younger children, p= ' num2str(p2) '.'])
else
    disp(['No difference in vertical tracts between adults and younger children, p= ' num2str(p2) '.'])
end


disp('=======================')
disp('Exploratory Contrasts')
disp(' ')
% Check for children (combining old and young) vs. adult.
% Re-assign groups.
gp_age_temp = gp_age;
gp_age_temp(find(gp_age_temp == 2)) = 1;
gp_age_temp(find(gp_age_temp == 3)) = 2;
% FOLLOW-UP: Is WM in h different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age_temp==1 & gp_hv==1)), toi_big(find(gp_age_temp==2 & gp_hv==1)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==1))) > nanmean(toi_big(find(gp_age_temp==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is greater in adults than in children (combining young and old), p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==1))) < nanmean(toi_big(find(gp_age_temp==1 & gp_hv==1)))
    disp([wm_measure ' in horizontal tracts is lower in adults than in children (combining young and old), p= ' num2str(p2) '.'])
else
    disp(['No difference in horizontal tracts between adults and children (combining young and old), p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in v different in children than in adults? no
[h, p2, ci] = ttest2(toi_big(find(gp_age_temp==1 & gp_hv==2)), toi_big(find(gp_age_temp==2 & gp_hv==2)));
if p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==2))) > nanmean(toi_big(find(gp_age_temp==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is greater in adults than in children (combining young and old), p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==2))) < nanmean(toi_big(find(gp_age_temp==1 & gp_hv==2)))
    disp([wm_measure ' in vertical tracts is lower in adults than in children (combining young and old), p= ' num2str(p2) '.'])
else
    disp(['No difference in in vertical tracts between adults and children (combining young and old), p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in children different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age_temp==1 & gp_hv==1)), toi_big(find(gp_age_temp==1 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==1 & gp_hv==1))) > nanmean(toi_big(find(gp_age_temp==1 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in children (combining young and old), p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==1 & gp_hv==1))) < nanmean(toi_big(find(gp_age_temp==1 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in children (combining young and old), p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in children (combining young and old), p= ' num2str(p2) '.'])
end
% FOLLOW-UP: Is WM in adults different in h than in v? yes
[h, p2, ci] = ttest2(toi_big(find(gp_age_temp==2 & gp_hv==1)), toi_big(find(gp_age_temp==2 & gp_hv==2))); % note that a two-sample ttest was used here because of unequal sample sizes
if p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==1))) > nanmean(toi_big(find(gp_age_temp==2 & gp_hv==2)))
    disp([wm_measure ' is greater in horizontal than vertical in adults, p= ' num2str(p2) '.'])
elseif p2 <= .05/mc && nanmean(toi_big(find(gp_age_temp==2 & gp_hv==1))) < nanmean(toi_big(find(gp_age_temp==2 & gp_hv==2)))
    disp([wm_measure ' is lower in horizontal than vertical in adults, p= ' num2str(p2) '.'])
else
    disp(['No difference in ' wm_measure ' between horizontal and vertical in adults, p= ' num2str(p2) '.'])
end
