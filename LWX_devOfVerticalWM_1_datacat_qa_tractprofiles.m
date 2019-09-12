% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

prefix = ''; % 'streamlinecount25920-'
blprojectid = ['/' prefix 'proj-5a74d1e26ed91402ce400cca/'];

w_measures = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'};

c = colorcube;
yc_color = [0 0.4470 0.7410 0.2];
oc_color = [0.4660 0.6740 0.1880 0.2];
a_color = [0.6350 0.0780 0.1840 0.2];

save_figures = 'yes';

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
addpath(genpath(fullfile(rootDir, 'proj-5a74d1e26ed91402ce400cca/')));

% Identify outliers to be removed.
% 128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
% 315 because strange WM in occipital lobe leading to no right or left VOF
% 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
outlier = [128 315 318]; %NOTE: need to update to not plot the outliers

fcount = 0;
for w = 1:length(w_measures)
    
    wm_measure = w_measures{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    sub_count = 0;
    for i = 1:size(grp_contents, 1)
        
        % Display current sub ID.
        disp(grp_contents(i).name)
        
        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-tractprofile*', 'profiles', '*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        for j = 1:size(sub_contents_tractprofiles)
            
            
            % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
            if i == 1 && j == 1
                
                tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                
            end
            
            % Read in data for this subject and this tract.
            data_temp = readtable(fullfile(sub_contents_tractprofiles(j).folder, sub_contents_tractprofiles(j).name));
            
            % Get middle 80%.
            start = size(data_temp, 1)*.1;
            stop = size(data_temp, 1)*.9;
            
            % Read in mean WM measure.
            if strcmp(wm_measure, 'ad')
                
                m_wm(:, j, sub_count) = data_temp.ad_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.ad_2(start:stop);
                
            elseif strcmp(wm_measure, 'fa')
                
                m_wm(:, j, sub_count) = data_temp.fa_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.fa_2(start:stop);
                
            elseif strcmp(wm_measure, 'md')
                
                m_wm(:, j, sub_count) = data_temp.md_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.md_2(start:stop);
                
            elseif strcmp(wm_measure, 'rd')
                
                m_wm(:, j, sub_count) = data_temp.rd_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.rd_2(start:stop);
                
            elseif strcmp(wm_measure, 'icvf')
                
                m_wm(:, j, sub_count) = data_temp.icvf_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.icvf_2(start:stop);
                
            elseif strcmp(wm_measure, 'isovf')
                
                m_wm(:, j, sub_count) = data_temp.isovf_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.isovf_2(start:stop);
                
            elseif strcmp(wm_measure, 'od')
                
                m_wm(:, j, sub_count) = data_temp.od_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.od_2(start:stop);
                
            end
            
            % Grab tract name for grouping variable.
            tract(:, j, sub_count) = repmat({sub_contents_tractprofiles(j).name(1:end-13)}, 161, 1);
            
            % Grab subID.
            sub(:, j, sub_count) = repmat(str2num(grp_contents(i).name(end-2:end)), 161, 1);
            
            clear data_temp
            
        end % sub_contents
        
    end % group_contents
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract));
    tract(t) = {'empty'};
    
    % Get group indices for tract names.
    G = findgroups(tract(:));
    for i = 1:length(sub)
        if sub(i) < 200
            group(i) = 1;
        elseif sub(i) < 300
            group(i) = 2;
        else
            group(i) = 3;
        end
    end
    group = group';
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get a list of unique sub IDs.
    subID = unique(sub);
    
    % Plot tract profiles for each tract.
    for k = 1:size(list_tract, 1)
        
        if ~strcmp(list_tract{k}, 'empty')
            
            % Find entries that are for this tract.
            t_idx = strcmp(tract, list_tract{k});
            
            % Open a new figure for this tract.
            fcount = fcount + 1;
            figure(fcount)
            
            yc_count = 0; oc_count = 0; a_count = 0;
            for s = 1:size(sub, 3)
                
                if subID(s) ~=0 && ~ismember(subID(s), outlier)
                    
                    % Find entries that are for this subject.
                    s_idx = sub == subID(s);
                    
                    % Subset the thing so that we only plot for this tract and for this subject.
                    t_temp = m_wm(find(t_idx == 1 & s_idx == 1));
                    
                    if ~isempty(t_temp)
                        
                        % Code the plot for subject.
                        if subID(s) < 200
                            
                            yc_count = yc_count + 1;
                            
                            % Young child.
                            plot(t_temp, 'LineStyle', '-', 'Color', yc_color)
                            
                            % Collect.
                            yc(:, yc_count) = t_temp;
                            
                        elseif subID(s) < 300
                            
                            oc_count = oc_count + 1;
                            
                            % Older child.
                            plot(t_temp, 'LineStyle', '-', 'Color', oc_color)
                            
                            % Collect.
                            oc(:, oc_count) = t_temp;
                            
                        else
                            
                            a_count = a_count + 1;
                            
                            % Adult.
                            plot(t_temp, 'LineStyle', '-', 'Color', a_color)
                            
                            % Collect.
                            a(:, a_count) = t_temp;
                            
                        end
                        hold on;
                        
                    end % if subID{s} ~= 0
                    
                end %if exist
                
            end %sub
            
            % Plot means and standard deviations.
            plot(mean(yc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', yc_color(1:3))
            hi = mean(yc, 2) + std(yc, 0, 2); lo = mean(yc, 2) - std(yc, 0, 2); x = (1:size(mean(yc, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], yc_color(1:3));
            set(hp1, 'facecolor', yc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(mean(oc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', oc_color(1:3))
            hi = mean(oc, 2) + std(oc, 0, 2); lo = mean(oc, 2) - std(oc, 0, 2); x = (1:size(mean(oc, 2),1))';
            hp2 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], oc_color(1:3));
            set(hp2, 'facecolor', oc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(mean(a, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', a_color(1:3))
            hi = mean(a, 2) + std(a, 0, 2); lo = mean(a, 2) - std(a, 0, 2); x = (1:size(mean(a, 2),1))';
            hp3 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], a_color(1:3));
            set(hp3, 'facecolor', a_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            ylabel(wm_measure);
            xlabel(list_tract{k});
            box off;
            
            print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-depsc')
            
            hold off;
            
        end % if
        
    end %tract
    
end % for w
    
    
    
    
    
    
    
    
