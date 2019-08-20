% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

% blprojectid = '/streamlinecount25920-proj-5a74d1e26ed91402ce400cca/';
blprojectid = '/proj-5a74d1e26ed91402ce400cca/';

b_measures = {'age', 'lit', 'vm', 'fm'};
w_measures = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'};

list_wm_toi = {'leftAslant', 'rightAslant', 'leftVOF', 'rightVOF', 'leftTPC', 'rightTPC', ...
    'leftMDLFspl', 'rightMDLFspl', 'leftMDLFang', 'rightMDLFang', 'leftpArc', 'rightpArc', ...
    'leftArc', 'rightArc', 'leftSLF1And2', 'rightSLF1And2', 'leftCST', 'rightCST'};

test_type = 'standard'; % standard, partial
pval = .05; p_bonf = pval/60; %length(list_wm_toi);
remove_outliers = 'yes';
save_figures = 'yes';

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
addpath(genpath([rootDir blprojectid]));

% Read in behavioral data.
load([rootDir 'supportFiles/LWX_all_groupings.mat']);

% Header for LWX Behavioral DatFa. gp assignments are
% young children (=1), older children (=2), and adults (=3).
% gp_lit (early literate = 1, literate = 2),
% gp_vm (low cm = 1, high vm = 2), and
% gp_fm (low cm = 1, high vm = 2).
header = {'subID', 'age_mo', 'vmi', 'vp', 'mc', 'pegs_dom', 'pegs_ndom', 'lwi', 'spell', 'wa', 'sos', 'c_vm', 'c_fm', 'c_lit', 'gp_age', 'gp_lit', 'gp_vm', 'gp_fm'};

fcount = 0;
for w = 1%:length(w_measures)
       
    for b = 1%:length(b_measures)
        
        wm_measure = w_measures{w};
        beh_measure = b_measures{b};
        
        %% BEHAVIOR.
        
        % Set names for correct columns given user input..
        if strcmp(beh_measure, 'age')
            
            measure_in = 'age_mo';
            gp_in = 'gp_age';
            
            outlier = [108 128 318];
            %128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
            % 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
            
        elseif strcmp(beh_measure, 'lit')
            
            measure_in = 'c_lit';
            gp_in = 'gp_lit';
            
            outlier = [108 128 318];
            
        elseif strcmp(beh_measure, 'vm')
            
            measure_in = 'c_vm';
            gp_in = 'gp_vm';
            
            outlier = [108 128 318];
            
        elseif strcmp(beh_measure, 'fm')
            
            measure_in = 'c_fm';
            gp_in = 'gp_fm';
            
            outlier = [108 128 318];
            
        end
        
        %% TRACTOGRAPHY.
        
        % Get contents of the directory where the tract measures for this subject are stored.
        grp_contents = dir([rootDir filesep blprojectid filesep]);
        
        % Remove the '.' and '..' files.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
        
        % Keep only names that are subject folders.
        grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
        
        % Load in each tract's tractography measures for this subject.
        sub_count = 0;
        for i = 1:size(grp_contents, 1)
            
            %     % Only read in data for subjects that we want to consider in this
            %     analysis. Right now we don't have all behavioral data entered.
            if sum(ismember(data_lwx(:, 1), grp_contents(i).name)) ~= 0
                
                % Display current sub ID.
                disp(grp_contents(i).name)
                
                % Update subject counter for when not all subjects are used/needed.
                sub_count = sub_count + 1;
                
                % Get contents of the directory where the tract measures for this subject are stored.
                sub_contents_tractprofiles = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-tractprofile*/profiles/*.csv']);
                sub_contents_tractstats = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-tractmeasures*/*.csv']);
                
                % Get contents of the directory where the SNR values for this subject are stored.
                sub_contents_snr = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-raw.tag-snr*/*product.json']);
                
                % Remove the '.' and '..' files.
                sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
                sub_contents_tractstats = sub_contents_tractstats(arrayfun(@(x) x.name(1), sub_contents_tractstats) ~= '.');
                sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
                
                % Read in stats for this subject's tracts.
                data_stats_temp = readtable([sub_contents_tractstats.folder filesep sub_contents_tractstats.name]);
                
                for j = 1:size(sub_contents_tractprofiles)
                    
                    % Read in data for this subject and this tract.
                    data_temp = csvread([sub_contents_tractprofiles(j).folder filesep sub_contents_tractprofiles(j).name], 1, 0);
                    
                    % Get middle 80%.
                    start = size(data_temp, 1)*.1;
                    stop = size(data_temp, 1)*.9;
                    
                    % Read in mean WM measure.
                    if strcmp(wm_measure, 'ad')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 1)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 2));
                        
                    elseif strcmp(wm_measure, 'fa')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 3)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 4));
                        
                    elseif strcmp(wm_measure, 'md')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 5)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 6));
                        
                    elseif strcmp(wm_measure, 'rd')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 7)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 8));
                        
                    elseif strcmp(wm_measure, 'icvf')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 17)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 18));
                        
                    elseif strcmp(wm_measure, 'isovf')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 19)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 20));
                        
                    elseif strcmp(wm_measure, 'od')
                        
                        m(sub_count, j) = mean(data_temp(start:stop, 21)); %sd_fa(sub_count, j) = mean(data_temp(start:stop, 22));
                        
                    end
                    
                    % Grab tract name for grouping variable.
                    tract{sub_count, j} = sub_contents_tractprofiles(j).name(1:end-13);
                    
                    % Grab subID.
                    sub(sub_count) = str2num(grp_contents(i).name(end-2:end));
                    
                    clear data_temp
                    
                    % Deal with differences in puncuation, spelling, etc.
                    if strcmp(tract{sub_count, j}, 'left_AnterioFrontoCerebellar')
                        
                        name_temp = strcat(tract{sub_count, j}(1:4), [' ' tract{sub_count, j}(6:end)]);
                        
                    elseif strcmp(tract{sub_count, j}, 'right_AnterioFrontoCerebellar')
                        
                        name_temp = strcat(tract{sub_count, j}(1:5), [' ' tract{sub_count, j}(7:end)]);
                        
                    elseif strncmp(tract{sub_count, j}, 'left_', 5)
                        
                        name_temp = strcat(tract{sub_count, j}(1:4), tract{sub_count, j}(6:end));
                        
                    elseif strncmp(tract{sub_count, j}, 'right_', 6)
                        
                        name_temp = strcat(tract{sub_count, j}(1:5), tract{sub_count, j}(7:end));
                        
                    else
                        
                        name_temp = tract{sub_count, j};
                        
                    end
                    
                    % Get number of streamlines for this tract and for this subject.
                    n_streamlines(sub_count, j) = data_stats_temp.StreamlineCount(find(strcmp(data_stats_temp.TractName, name_temp)));
                    
                    clear name_temp
                    
                end % end j
                
            end % ismember
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            for g = 1:size(data_snr_temp.SNRInB0_X_Y_Z, 1)
                
                get_temp(g) = str2num(data_snr_temp.SNRInB0_X_Y_Z{g});
                
            end
            
            snr(sub_count) = min(get_temp);
            
            clear data_stats_temp data_snr_temp get_temp
            
        end
        
        % Find empty cells.
        t = find(cellfun(@isempty,tract));
        
        % Enter 'empty' in empty cells.
        tract(t) = {'empty'};
        
        % Get indexing for tract groups based on tract names.
        G = findgroups(tract(:));
        
        % If this is the first time through this loop, plot the basic FA
        % scatter plots with no behavioral correlate and bar plot for number of streamlines per tract, for QA.
        if b==1 && w==1
            
            if strcmp(wm_measure, 'fa')
                
                % Visualize: fa.
                fcount = fcount + 1;
                figure(fcount)
                x = repmat(transpose(1:size(m, 1)), [1 size(m, 2)]);
                colors = colorcube(size(m, 2));
                gscatter(x(:), m(:), G, colors)
                b_plot = gca; legend(b_plot,'off');
                
                ylabel('Mean FA per Tract (averaged across streamlines)')
                xlabel('SubID')
                
                print([rootDir 'plots/plot_scatter_' wm_measure], '-dpng')
                print([rootDir 'plots/eps/plot_scatter_' wm_measure], '-depsc')
                
            end
            
            % Visualize: n_streamlines.
            fcount = fcount + 1;
            figure(fcount)
            x = repmat(transpose(1:size(n_streamlines, 1)), [1 size(n_streamlines, 2)]);
            colors = colorcube(size(n_streamlines, 2));
            gscatter(x(:), n_streamlines(:), G, colors)
            b_plot = gca; legend(b_plot,'off');
            
            ylabel('Number of Streamlines per Tract')
            xlabel('SubID')
            
            print([rootDir 'plots/plot_scatter_n_streamlines'], '-dpng')
            print([rootDir 'plots/eps/plot_scatter_n_streamlines'], '-depsc')
            
            fcount = fcount + 1;
            figure(fcount)
            bar(nanmean(n_streamlines, 1))
            ylabel('Mean Number of Streamlines per Tract (across subjects)')
            xlabel('Tract')
            
            print([rootDir 'plots/plot_bar_n_streamlines'], '-dpng')
            print([rootDir 'plots/eps/plot_bar_n_streamlines'], '-depsc')
            
        end
        
        %% ================== Relate WM microstructure to age and behavioral measures ================== %%
        
        % Set the mean of any tracts with less than 100 streamlines to NaN so that they will not be included in the average.
        m(n_streamlines <= 100) = NaN;
        
        % Get a list of unique tract names.
        list_tract = unique(tract);
        
        % Get WM measurements for each tract (reorganizing so that each column is a tract).
        for k = 1:size(list_tract, 1)
            
            % Get the mean of the wm_measure for this tract (take sum because don't want to include zero columns; only one value of interest per row).
            y_temp(:, k) = nansum(m.*strcmp(tract, list_tract{k}), 2);
            
        end % end k
        y = cat(2, transpose(sub), y_temp); clear y_temp
        
        % Create new grouping vectors, one subject at a time.
        for s = 1:length(sub)
            
            if ismember(sub(s), data_lwx(:, find(strcmp(header, 'subID'))))
                
                % Get this subject's AGE group.
                group(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, gp_in)));
                
                % Get this subject's MEASURE.
                measure(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, measure_in)));
                
                % Get this subject's AGE.
                cov_age(s) = data_lwx(find(sub(s) == data_lwx(:, find(strcmp(header, 'subID')))), find(strcmp(header, 'age_mo')));
                
                % Get this subject's WM measurements.
                wm(s, :) = y(find(y(:, 1) == sub(s)), 3:end); % start at 3 because first column is sub and second column is 'empty' (i.e., 0)
                
            end
            
        end % end s
        
        % Use only subjects for whom we have both behavioral and tractography measures.
        
        % Remove the 'group' of zero (i.e., empty). Make sure to modify 'group' last.
        measure = transpose(measure(group~=0));
        n_streamlines = n_streamlines(group~=0, :);
        wm = wm(group~=0, :);
        cov_age = transpose(cov_age(group~=0));
        sub = transpose(sub(group~=0));
        snr = transpose(snr(group~=0));
        group = transpose(group(group~=0));
        
        % Set all zeros to NaN.
        measure(measure == 0) = NaN; wm(wm == 0) = NaN; cov_age(cov_age == 0) = NaN; sub(sub == 0) = NaN; n_streamlines(n_streamlines == 0) = NaN; snr(snr == 0) = NaN; group(group == 0) = NaN;
        
        % Remove outliers.
        if strcmp(remove_outliers, 'yes') && exist('outlier')
            
            idx_remove = ~any(~(sub~=outlier), 2);
            
            wm = wm(idx_remove, :);
            n_streamlines = n_streamlines(idx_remove, :);
            measure = measure(idx_remove);
            cov_age = cov_age(idx_remove);
            group = group(idx_remove);
            snr = snr(idx_remove);
            sub = sub(idx_remove);
            
        end
        
        % Make a data set for only the child data. 
        wm_childrenOnly = wm(group~=3, :);
        n_streamlines_childrenOnly = n_streamlines(group~=3, :);
        measure_childrenOnly = measure(group~=3);
        cov_age_childrenOnly = cov_age(group~=3);
        sub_childrenOnly = sub(group~=3);
        snr_childrenOnly = snr(group~=3);
        group_childrenOnly = group(group~=3);
        
        % Are there AGE group differences in FA?
        for test = 1:size(wm, 2)
              
            % One-way ANOVA, display turned off
            [p, tbl, stats] = anova1(wm_childrenOnly(:, test), group_childrenOnly, 'off');
                
            % Output information to the command window.
            if p <= p_bonf
                
                disp([list_tract{test} '     ' num2str(p) '     *, passed MC'])
                
            elseif p <= pval
                
                disp([list_tract{test} '     ' num2str(p) '     *'])
                
            elseif p <= 0.1
                
                disp([list_tract{test} '     ' num2str(p)])
                
            end
            
        end % end test
        
        disp('==========================')
        
        % Z-score measure for correlations, if not age.
        if ~strcmp(measure_in, 'age_mo')
            
            idxy = isnan(measure).*isnan(measure);
            measure_z = (measure - nanmean(measure, 1))./nanstd(measure, 1);
            measure_childrenOnly_z = (measure_childrenOnly - nanmean(measure_childrenOnly, 1))./nanstd(measure_childrenOnly, 1);
            
        else
            
            measure_z = measure;
            measure_childrenOnly_z = measure_childrenOnly;
            
        end
        
        % Zscore WM based on within-tract means for correlations
        idxy = isnan(wm).*isnan(wm);
%         wm_z = (wm - nanmean(wm, 1))./nanstd(wm, 1);
%         wm_childrenOnly_z = (wm_childrenOnly - nanmean(wm_childrenOnly, 1))./nanstd(wm_childrenOnly, 1);
        wm_z = (wm - nanmean(wm(:)))./nanstd(wm(:));
        wm_childrenOnly_z = (wm_childrenOnly - nanmean(wm_childrenOnly(:)))./nanstd(wm_childrenOnly(:));
        
        % SECOND: Does FA correlate with age in children?
        for test = 1:size(wm, 2)
            
            if strcmp(test_type, 'standard')
                
                [r, p] = corrcoef(cat(2, measure_childrenOnly_z, wm_childrenOnly_z(:, test)), 'Rows', 'pairwise');
                p_check = p(2, 1);
                
            elseif strcmp(test_type, 'partial')
                
                [r, p] = partialcorr(measure_childrenOnly_z, wm_childrenOnly_z(:, test), cov_age_childrenOnly, 'Rows', 'pairwise');
                p_check = p;
                
            end
            
            % Output information to the command window.
            if p_check <= p_bonf
                
                disp([list_tract{test} '     ' num2str(p_check) '     *, passed MC'])
                
            elseif p_check <= pval
                
                disp([list_tract{test} '     ' num2str(p_check) '     *'])
                
            elseif p_check <= 0.1
                
                disp([list_tract{test} '     ' num2str(p_check)])
                
            end
            
            % Only plot if passes minimal stastical threshold and if it's a standard correlation.
            %     if p_check <= 0.1 && strcmp(test_type, 'standard')
            
            % Only plot if this is a tract of interest (toi).
            if sum(contains(list_wm_toi, list_tract{test})) > 0
                
                fcount = fcount + 1;
                
                figure(fcount)
                
                plot(measure_childrenOnly_z(group_childrenOnly==1), wm_childrenOnly_z(group_childrenOnly==1, test), 'bo')
                hold on
                plot(measure_childrenOnly_z(group_childrenOnly==2), wm_childrenOnly_z(group_childrenOnly==2, test), 'go')
                
                if strcmp(beh_measure, 'age')
                    plotcorr(measure_childrenOnly_z, wm_childrenOnly_z(:, test), [beh_measure ' (mo)'], [list_tract{test} ', average ' wm_measure])
                else
                    plotcorr(measure_childrenOnly_z, wm_childrenOnly_z(:, test), beh_measure, [list_tract{test} ', average ' wm_measure])
                end
                
                % Identify Minimum FA subject.
                idy_minFA = min(wm_childrenOnly_z(:, test));
                idx_minFA = measure_childrenOnly_z(find(wm_childrenOnly_z(:, test) == idy_minFA));
                sub_minFA = sub_childrenOnly(find(wm_childrenOnly_z(:, test) == idy_minFA));
                text(idx_minFA, idy_minFA, ['Sub ' num2str(sub_minFA) '(' num2str(idx_minFA) ',' num2str(idy_minFA) ')'])
                
                % Identify Maximum FA subject.
                idy_maxFA = max(wm_childrenOnly_z(:, test));
                idx_maxFA = measure_childrenOnly_z(find(wm_childrenOnly_z(:, test) == idy_maxFA));
                sub_maxFA = sub_childrenOnly(find(wm_childrenOnly_z(:, test) == idy_maxFA));
                text(idx_maxFA, idy_maxFA, ['Sub '  num2str(sub_maxFA) '(' num2str(idx_maxFA) ',' num2str(idy_maxFA) ')'])
                
                % Save plots, if requested.
                if strcmp(save_figures, 'yes') || strcmp(remove_outliers, 'yes')
                    
                    if strcmp(beh_measure, 'age')
                        ylim([-3, 3]); xlim([50, 110]);
                    end
                    
                    print([rootDir 'plots/plot_scatter_' beh_measure '_' wm_measure '_' list_tract{test} '_tractz'], '-dpng')
                    print([rootDir 'plots/eps/plot_scatter_' beh_measure '_' wm_measure '_' list_tract{test} '_tractz'], '-depsc')
                    
                elseif strcmp(save_figures, 'yes')
                    
                    print([rootDir 'plots/plot_scatter_' beh_measure '_' wm_measure '_' list_tract{test} '_wOutliers'], '-dpng')
                    print([rootDir 'plots/eps/plot_scatter_' beh_measure '_' wm_measure '_' list_tract{test} '_wOutliers'], '-depsc')
                    
                end
                
                hold off;
                
                % Plot Number of Streamlines per Subject
                
                fcount = fcount + 1;
                
                figure(fcount)
                
                plot(measure_childrenOnly_z(group_childrenOnly==1), n_streamlines_childrenOnly(group_childrenOnly==1, test), 'bo')
                hold on
                plot(measure_childrenOnly_z(group_childrenOnly==2), n_streamlines_childrenOnly(group_childrenOnly==2, test), 'go')
                
                if strcmp(beh_measure, 'age')
                    plotcorr(measure_childrenOnly_z, n_streamlines_childrenOnly(:, test), [beh_measure ' (mo)'], [list_tract{test} ', N streamlines'])
                else
                    plotcorr(measure_childrenOnly_z, n_streamlines_childrenOnly(:, test), beh_measure, [list_tract{test} ', N streamlines'])
                end
                
                % Identify Minimum FA subject.
                idy_minFA = min(n_streamlines_childrenOnly(:, test));
                idx_minFA = measure_childrenOnly_z(find(n_streamlines_childrenOnly(:, test) == idy_minFA));
                sub_minFA = sub_childrenOnly(find(n_streamlines_childrenOnly(:, test) == idy_minFA));
                text(idx_minFA, idy_minFA, ['Sub ' num2str(sub_minFA) '(' num2str(idx_minFA) ',' num2str(idy_minFA) ')'])
                
                % Identify Maximum FA subject.
                idy_maxFA = max(n_streamlines_childrenOnly(:, test));
                idx_maxFA = measure_childrenOnly_z(find(n_streamlines_childrenOnly(:, test) == idy_maxFA));
                sub_maxFA = sub_childrenOnly(find(n_streamlines_childrenOnly(:, test) == idy_maxFA));
                text(idx_maxFA, idy_maxFA, ['Sub ' num2str(sub_maxFA) '(' num2str(idx_maxFA) ',' num2str(idy_maxFA) ')'])
                
                % Save plots, if requested.
                if strcmp(save_figures, 'yes') || strcmp(remove_outliers, 'yes')
                    
                    print([rootDir 'plots/plot_scatter_' beh_measure '_nstreamlines_' list_tract{test} '_tractz'], '-dpng')
                    print([rootDir 'plots/eps/plot_scatter_' beh_measure '_nstreamlines_' list_tract{test} '_tractz'], '-depsc')
                    
                elseif strcmp(save_figures, 'yes')
                    
                    print([rootDir 'plots/plot_scatter_' beh_measure '_nstreamlines_' list_tract{test} '_wOutliers'], '-dpng')
                    print([rootDir 'plots/eps/plot_scatter_' beh_measure '_nstreamlines_' list_tract{test} '_wOutliers'], '-depsc')
                    
                end
                
                hold off
                
            end
            
        end % end test
        
        % Save all variables.
        save([rootDir 'supportFiles/LWX_data_' wm_measure '_' beh_measure '_tractz.mat'])
        
        % If this is the first time through, plot snr, for qa.
        if b==1 && w==1
            
            fcount = fcount + 1;
            figure(fcount)
            plot(measure_z(group==1), snr(group==1), 'bo') % use all subjects for this
            hold on
            plot(measure_z(group==2), snr(group==2), 'go')
            plot(measure_z(group==3), snr(group==3), 'ro')
            
            if strcmp(beh_measure, 'age')
                plotcorr(measure_z(group~=3), snr(group~=3), [beh_measure ' (mo)'], 'SNR')
            else
                plotcorr(measure_z, snr, beh_measure, 'SNR')
            end
            
            print([rootDir 'plots/plot_scatter_snr'], '-dpng')
            print([rootDir 'plots/eps/plot_scatter_snr'], '-depsc')
            
        end
        
        % Reset for next loop.      
        clearvars -except w b rootDir header fcount blprojectid b_measures w_measures list_wm_toi test_type pval p_bonf remove_outliers save_figures data_lwx
        
    end
    
end






