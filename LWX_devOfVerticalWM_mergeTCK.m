% This script will read in all available TCK files for a subject and merge
% them into one TCK file. BIDS format is assumed. 

% Dependencies: 
% mrtrix3

clear all; close all; clc;

rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
blprojectid = 'proj-5a74d1e26ed91402ce400cca';

% % read in tck file
% test_tck = read_mrtrix_tracks('/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/mergedtcks/sub-103_track.tck');

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
for i = 1:size(grp_contents, 1)
    
    % Display current sub ID.
    disp(grp_contents(i).name)
    
    % Get contents of the directories where the .tck files for this subject are stored.
    sub_contents_tck = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-track-tck.id*'));
    
    % Remove the '.' and '..' files.
    sub_contents_tck = sub_contents_tck(arrayfun(@(x) x.name(1), sub_contents_tck) ~= '.');
    
    % Get subID.
    subID = sub_contents_tck(1).folder(strfind(sub_contents_tck(1).folder, 'sub'):end);
    
    % Read in all the .tck files available for each subject, one at a time.
    for j = 1:length(sub_contents_tck)
        
        % Read in data for this subject and this tck.
        tck_temp = read_mrtrix_tracks(fullfile(sub_contents_tck(j).folder, sub_contents_tck(j).name, 'track.tck'));
        
        % Merge.
        if j == 1
            
            % Initiate new tck information for this subject.
            mergedtck = tck_temp.data;
            
            % Initiate seed and track counters for this subject.
            count_seed = tract_temp.max_num_seeds;
            count_track = tract_temp.max_num_tracks;
            
        else
            
            % Append tck for this tck to previously added tcks for this subject.
            mergedtck = [mergedtck tck_temp.data];
            
            % Update counters for seed and track counts for this subject.
            count_seed = count_seed + tract_temp.max_num_seeds;
            count_track = count_track + tract_temp.max_num_tracks;
            
        end
        
    end
    
    % Replace the streamline data in the track_temp struct so that we keep the header information. 
    tck_temp.data = mergedtck;
    
    % Update basic header information -- NOTE: Does this result in accurate header information for the merged tck?
    tck_temp.count = length(tck_temp.data);
    tck_temp.total_count = length(tck_temp.data);
    
    % Write out merged tck file.
    write_mrtrix_tracks(tck_temp, fullfile(rootDir, 'mergedtcks', [subID '_track.tck']))
    
    % Restart. 
    clear tck_temp mergedtck count_seed count_track sub_contents_tck subID
    
end


