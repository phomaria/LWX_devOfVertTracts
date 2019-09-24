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

wm_measure_here = {'fa', 'ad', 'md', 'rd', 'od'}; %, 'icvf', 'isovf'
% lo = [.49 1.37 .81 .52 .21 .65 .11];
% mid = [.53 1.43 .87 .59 .23 .75 .13];
% hi = [.57 1.49 .93 .66 .25 .85 .15];
lo = [.45 1.05 .65 .44 .12];
mid = [.50 1.15 .75 .52 .13];
hi = [.55 1.25 .85 .60 .14];

capsize = 0;
marker = 'o';
markeredgecolor_h = [0 .73 .73]'; markeredgecolor_v = [.146 0 0]';
markerfacecolor_h = [0 .73 .73]; markerfacecolor_v = [.146 0 0];
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0.05;
xticklength = 0;
xtickvalues = [1 2 3];
xlim_lo = 0.5; xlim_hi = 3.5;
save_figures = 'yes';

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    data_tbl = readtable(fullfile(rootDir, 'supportFiles', ['LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '_singleshell.csv']));
    
    % Convert into array and header for ease.
    data_all_in_header = data_tbl.Properties.VariableNames;
    data_all_in = table2array(data_tbl);

    % Get indices of subjects whose white matter values are all NaN.
    idx_notnan = ~all(isnan(data_all_in(:, 8:end)), 2);
    
    % Remove subjects whose white matter values are all NaN from the table and from the array.
    data_tbl = data_tbl(idx_notnan, :);
    data_all_in = data_all_in(idx_notnan, :);
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Get easy index for age group.
    group = data_tbl.group_age;

    % Plot means. %NOTE: Change errorbar to plot so that I can do whiskers instead of T-bars.
    figure(w)
    errorbars = errorbar(xtickvalues, [mean(data_tbl.meanH(group == 1)) mean(data_tbl.meanH(group == 2)) mean(data_tbl.meanH(group == 3))], ...
        [std(data_tbl.meanH(group == 1), 0, 1)/sqrt(length(data_tbl.meanH(group == 1))), std(data_tbl.meanH(group == 2), 0, 1)/sqrt(length(data_tbl.meanH(group == 2))), std(data_tbl.meanH(group == 3), 0, 1)/sqrt(length(data_tbl.meanH(group == 3)))]);
    set(errorbars, 'Color', markeredgecolor_h, 'Marker', marker, 'MarkerEdgeColor', markeredgecolor_h, 'MarkerFaceColor', markerfacecolor_h, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    hold on;    
    errorbars = errorbar(xtickvalues, [mean(data_tbl.meanV(group == 1)) mean(data_tbl.meanV(group == 2)) mean(data_tbl.meanV(group == 3))], ...
        [std(data_tbl.meanV(group == 1), 0, 1)/sqrt(length(data_tbl.meanV(group == 1))), std(data_tbl.meanV(group == 2)/sqrt(length(data_tbl.meanV(group == 2))), 0, 1), std(data_tbl.meanV(group == 3), 0, 1)/sqrt(length(data_tbl.meanV(group == 3)))]);
    set(errorbars, 'Color', markeredgecolor_v, 'Marker', marker, 'MarkerEdgeColor', markeredgecolor_v, 'MarkerFaceColor', markerfacecolor_v, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
   
    % xaxis
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = xtickvalues;
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    xlabels = {'Younger,Children', '  Older,Children', 'Adults'};
    xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
    xax.TickLabels = xlabels;
    xax.FontName = fontname;
    xax.FontSize = fontsize;
    xax.FontAngle = fontangle;

    % yaxis
    yax = get(gca,'yaxis');
    yax.Limits = [lo(w) hi(w)];
    yax.TickValues = [lo(w) mid(w) hi(w)];
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    yax.TickLabels = {num2str(lo(w), '%1.2f'), '', num2str(hi(w), '%1.2f')};
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    
    % general
    a = gca;
%     a.TitleFontWeight = 'normal';
    box off
    if strcmp(wm_measure_here{w}, 'od')
        a.YLabel.String = 'Orientation Dispersion (OD)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northwest', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'icvf')
        a.YLabel.String = 'Neurite Density (ICVF)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'isovf')
        a.YLabel.String = 'Isotropic Volume Fraction (ISOVF)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'fa')
        a.YLabel.String = 'Fractional Anisotropy (FA)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'ad')
        a.YLabel.String = 'Axial Diffusivity (AD)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'md')
        a.YLabel.String = 'Mean Diffusivity (MD)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    elseif strcmp(wm_measure_here{w}, 'rd')
        a.YLabel.String = 'Radial Diffusivity (RD)';
        legend({'Horizontal', 'Vertical'}, 'Location', 'northeast', 'Orientation', 'vertical');
        legend('boxoff');
    end
    a.YLabel.FontSize = fontsize;
    pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

    % Write.
    if strcmp(save_figures, 'yes')
        
        print(fullfile(rootDir, 'plots', ['plot_anova_' wm_measure_here{w} '_hv_yoa_singleshell']), '-dpng')
        print(fullfile(rootDir, 'plots', 'eps', ['plot_anova_' wm_measure_here{w} '_hv_yoa_singleshell']), '-depsc')
        
    end
    
    hold off;
    
end


