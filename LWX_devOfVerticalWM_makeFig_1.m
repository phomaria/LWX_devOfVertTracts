% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';

wm_measure_here = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'}; 
lo = [.49 1.37 .81 .52 .21 .65 .11];
mid = [.53 1.43 .87 .59 .23 .75 .13];
hi = [.57 1.49 .93 .66 .25 .85 .15];

marker = 'o';
markeredgecolor_h = [0 .73 .73]'; markeredgecolor_v = [.146 0 0]';
markerfacecolor_h = 'none'; markerfacecolor_v = 'none';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
save_figures = 'yes';

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    data_tbl = readtable(fullfile(rootDir, 'supportFiles', ['LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '.csv']));
    
    % Convert into array and header for ease.
    data_all_in = table2array(data_tbl);
    data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Get easy index for age group.
    group = data_tbl.group_age;

    % Plot means.
    figure(w)
    errorbars = errorbar([1, 2, 3], [mean(data_tbl.meanH(group == 1)) mean(data_tbl.meanH(group == 2)) mean(data_tbl.meanH(group == 3))], ...
        [std(data_tbl.meanH(group == 1), 0, 1)/sqrt(length(data_tbl.meanH(group == 1))), std(data_tbl.meanH(group == 2), 0, 1)/sqrt(length(data_tbl.meanH(group == 2))), std(data_tbl.meanH(group == 3), 0, 1)/sqrt(length(data_tbl.meanH(group == 3)))]);
    set(errorbars,  'Color', markeredgecolor_h, 'Marker', marker, 'MarkerEdgeColor', markeredgecolor_h, 'MarkerFaceColor', markerfacecolor_h, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle);
    hold on;    
    errorbars = errorbar([1, 2, 3], [mean(data_tbl.meanV(group == 1)) mean(data_tbl.meanV(group == 2)) mean(data_tbl.meanV(group == 3))], ...
        [std(data_tbl.meanV(group == 1), 0, 1)/sqrt(length(data_tbl.meanV(group == 1))), std(data_tbl.meanV(group == 2)/sqrt(length(data_tbl.meanV(group == 2))), 0, 1), std(data_tbl.meanV(group == 3), 0, 1)/sqrt(length(data_tbl.meanV(group == 3)))]);
    set(errorbars,  'Color', markeredgecolor_v, 'Marker', marker, 'MarkerEdgeColor', markeredgecolor_v, 'MarkerFaceColor', markerfacecolor_v, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle);
    set(gca, 'xtick', [1 2 3]); set(gca, 'ytick', [lo(w) mid(w) hi(w)]);
    a = gca;
    xlabels = {'Younger Children', 'Older Children','Adults'};
    xlabels = cellfun(@(x) strrep(x,' ','\newline'), xlabels,'UniformOutput',false);
    a.XTickLabel = xlabels;
    a.YTickLabel = {num2str(lo(w), '%1.2f'), num2str(mid(w), '%1.2f'), num2str(hi(w), '%1.2f')};
    set(gca, 'TickLength', [0 0])
    set(gca, 'FontSize', fontsize);
    set(gca, 'FontName', fontname);
    xlim([0.5 3.5]); ylim([lo(w) hi(w)]);
    box off
%     legend({'Horizontal', 'Vertical'}, 'Location', 'northinside', 'Orientation', 'horizontal');
    
    if strcmp(wm_measure_here{w}, 'od')
        ylabel('Orientation Dispersion (OD)');
    elseif strcmp(wm_measure_here{w}, 'icvf')
        ylabel('Neurite Density (ICVF)');
    elseif strcmp(wm_measure_here{w}, 'isovf')
        ylabel('Isotropic Volume Fraction (ISOVF)');
    elseif strcmp(wm_measure_here{w}, 'fa')
        ylabel('Fractional Anisotropy (FA)');
    elseif strcmp(wm_measure_here{w}, 'ad')
        ylabel('Axial Diffusivity (AD)');
    elseif strcmp(wm_measure_here{w}, 'md')
        ylabel('Mean Diffusivity (MD)');
    elseif strcmp(wm_measure_here{w}, 'rd')
        ylabel('Radial Diffusivity (RD)');
    end

    if strcmp(save_figures, 'yes')
        
        print(fullfile(rootDir, 'plots', ['plot_anova_' wm_measure_here{w} '_hv_yoa']), '-dpng')
        print(fullfile(rootDir, 'plots', 'eps', ['plot_anova_' wm_measure_here{w} '_hv_yoa']), '-depsc')
        
    end
    
    hold off;
    
end


