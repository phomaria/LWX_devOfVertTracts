% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

remove_outliers = 'no';

prefix = ''; % 'streamlinecount25920-'
blprojectid = ['/' prefix 'proj-5a74d1e26ed91402ce400cca/'];

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Identify outliers to be removed.
% 128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
% 315 because strange WM in occipital lobe leading to no right or left VOF
% 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
outlier = [108];% 128 315 318];

%% WHITE MATTER MEASURES

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
load(fullfile(rootDir, 'supportFiles', 'LWX_data_fa_raw_singleshell.mat'))

% Convert into array and header for ease.
data_all_in = table2array(data_tbl);
data_all_in_header = data_tbl.Properties.VariableNames;

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
    
    % Get contents of the directory where the SNR values for this subject are stored.
    sub_contents_snr = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, '/dt-raw.tag-snr*/*product.json'));
    % Remove the '.' and '..' files.
    sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
    
    % Get SNR for this subject.
    data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
    
    for g = 1:size(data_snr_temp.SNRInB0_X_Y_Z, 1)
        
        get_temp(g) = str2num(data_snr_temp.SNRInB0_X_Y_Z{g});
        
    end
    
    % Get subID.
    subID(sub_count) = str2num(grp_contents(i).name(5:7));
    
    % Get SNR.
    snr(sub_count) = min(get_temp);
    
    % Get age group.
    group(sub_count) = str2num(grp_contents(i).name(5));
    
    % Get age in months.
    age(sub_count) = beh_data_in_tbl.Age_months(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(i).name(5:7)))));
    
    clear data_snr_temp get_temp
    
end

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    snr = snr(~idx_outlier);
    group = group(~idx_outlier);
    age = age(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
        % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

figure
plot(age(group==1)', snr(group==1)', 'bo') % use all subjects for this
hold on
plot(age(group==2)', snr(group==2)', 'go')
% plot(age(group==3)', snr(group==3)', 'ro')

plotcorr(age(group~=3)', snr(group~=3)', 'Age (mo)', 'SNR', ttlstr)

print([rootDir 'plots/plot_scatter_snr_singleshell'], '-dpng')
print([rootDir 'plots/eps/plot_scatter_snr_singleshell'], '-depsc')

hold off;

% Write out table for anova.
t_out = array2table(cat(2, subID', group', age', snr'), 'VariableNames', {'subID', 'group_age', 'cov_age', 'snr'});
writetable(t_out, fullfile(rootDir, 'supportFiles', 'LWX_devOfVerticalWM_forSPSS_SNR_singleshell.csv'));
    



