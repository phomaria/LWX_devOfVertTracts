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
wm_measure = 'fa'; %fa, od, icvf

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_tractz.mat'])
clearvars -except save_figures wm group list_tract wm_measure rootDir covariates sub beh_measure beh_measures

hv_list = {'leftSLF1And2' 'rightSLF1And2', 'leftIFOF', 'rightIFOF', 'leftILF', 'rightILF', 'leftArc', 'rightArc', 'leftSLF3', 'rightSLF3', ...
    'leftAslant', 'rightAslant', 'leftTPC', 'rightTPC', 'leftpArc', 'rightpArc', 'leftMDLFspl', 'rightMDLFspl', 'leftVOF', 'rightVOF', 'leftMDLFang', 'rightMDLFang'};

count = 0;
for s = 1:length(hv_list)
    
    % Get index matrices for hypothesis-driven grouping of WM tracts. LH only right now.
    for k = 1:length(list_tract)
        
        % Indices of horizontal tracts.
        hv(k) = strcmp(list_tract{k}, hv_list{s});
        
    end
    
    if sum(hv) ~= 0
        
        count = count + 1;
        
        % Select the measurements of the tracts that I care about and convert all zeros to NaN.
        toi = wm(:, find(hv ~= 0)); toi(toi==0) = NaN;
        
        % Update the tract indexing to correspond with toi dimensions.
        hv = hv(find(hv~=0));
        
        % Update the tract list so that we can grab the tract name.
        list_tract2 = list_tract(hv~=0);
        
        % Reorganize data to prepare for two-way repeated measures anova using anovan: dependent variable.
        toi_big = toi(:); %not needed
        
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
        % variable and then for the subject grouping variable (gp_age) to be noted as being nested
        % under the subID.
        % Test the 2-ways interactions.
        % [p, tbl, stats] = anovan(toi_big, {gp_age, gp_hv, gp_sub}, 'random', 3, 'nested', [0 0 0 ; 0 0 0 ; 1 0 0], ...
        %     'model','interaction', 'varnames',{'Age Group', 'Tract Orientation', 'subID'}, 'display', 'off');
        [p, tbl, stats] = anova1(toi_big, gp_age, 'off');
        
        % Looking for an interaction between Age Group and Tract Orientation,
        % so report on p(4), in this case.
        disp(p)
        
        % Visualize.
        figure(count)
        c = colorcube;
        yc_color = c(1:11, :);
        oc_color = c(12:23, :);
        a_color = c(24:35, :);
        G = findgroups(gp_sub);
        gradient = transpose(linspace(.5, 1, length(unique(G(find(gp_age==1))))));
        % Visualize: young children, horizontal
        h1 = gscatter(ones(size(find(gp_age==1), 1), 1), toi_big(find(gp_age==1)), G(find(gp_age==1)), yc_color);
        %'.', 'MarkerFaceColor', yc_color, 'MarkerEdgeColor', yc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
        hold on;
        plot(linspace(.8, 1.2, 5), repmat(nanmean(toi_big(find(gp_age==1 & gp_hv==1))), [1 5]), 'k-')
        
        % Visualize: older children, horizontal
        h2 = gscatter(2*ones(size(find(gp_age==2), 1), 1), toi_big(find(gp_age==2)), G(find(gp_age==2)), oc_color);
        %'.', 'MarkerFaceColor', oc_color, 'MarkerEdgeColor', oc_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
        plot(linspace(1.8, 2.2, 5), repmat(nanmean(toi_big(find(gp_age==2 & gp_hv==1))), [1 5]), 'k-')
        
        % Visualize: adults, horizontal
        h3 = gscatter(3*ones(size(find(gp_age==3), 1), 1), toi_big(find(gp_age==3)), G(find(gp_age==3)), a_color);
        %'.', 'MarkerFaceColor', a_color, 'MarkerEdgeColor', a_color, 'MarkerFaceAlpha', 6/8, 'MarkerEdgeAlpha', 6/8);
        plot(linspace(2.8, 3.2, 5), repmat(nanmean(toi_big(find(gp_age==3 & gp_hv==1))), [1 5]), 'k-')
        
        b_plot = gca; legend(b_plot,'off');
        
        % Dividing lines.
        buffer = 0.3; linesize = 2; fontsize= 10;
        
        set(gca,'xtick',[1 2 3]);
        % set(gca,'xticklabel',{'children', 'children', 'adults', 'children', 'children', 'adults'});...
        %     '4.5-6.0 yrs.', '6.0-8.5 yrs.', ' ', '4.5-6.0 yrs.', '6.0-8.5 yrs.', ' '});
        labels = {'Children 4.5-6.0yrs.','Children 6.0-8.5yrs.','Adults'};
        labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
        a = gca;
        a.XTickLabel = labels;
        set(gca,'FontSize', fontsize);
        
        if strcmp(wm_measure, 'fa')
            ylabel('Fractional Anisotropy');
            ylim([0.3 0.75]);
        elseif strcmp(wm_measure, 'od')
            ylabel('Orientation Dispersion');
            ylim([0.5 1.2]);
        elseif strcmp(wm_measure, 'icvf')
            ylabel('Neurite Density (ICVF)');
            ylim([0.2 0.9]);
        end
        xlim([0.5 3.5]);
        title([hv_list{s} ', F=' num2str(tbl{2, 5}) ', p=' num2str(tbl{2, 6})])
        
        box off;
        
        print([rootDir 'plots/plot_3x2anova_' wm_measure '_hv_' hv_list{s}], '-dpng')
        print([rootDir 'plots/eps/plot_3x2anova_' wm_measure '_hv_' hv_list{s}], '-depsc')
        
        hold off;
        
    end
    
end


