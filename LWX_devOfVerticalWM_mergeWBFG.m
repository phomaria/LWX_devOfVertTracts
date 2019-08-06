clear all; close all; clc;

rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';
blprojectid = 'proj-5a74d1e26ed91402ce400cca';
%
% % read in tck file
% track_test = read_mrtrix_tracks('/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/mergedtcks/sub-103_track.tck');

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir([rootDir blprojectid filesep]);

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
for i = 1:size(grp_contents, 1)
    
    % Display current sub ID.
%     disp(grp_contents(i).name)
    
    % Get contents of the directories where the .tck files for this subject are stored.
    sub_contents_tck = dir([grp_contents(i).folder filesep grp_contents(i).name '/dt-neuro-track-tck.id*/']);
    
    % Remove the '.' and '..' files.
    sub_contents_tck = sub_contents_tck(arrayfun(@(x) x.name(1), sub_contents_tck) ~= '.');
    
    % Get subID.
    subID = sub_contents_tck(1).folder(strfind(sub_contents_tck(1).folder, 'sub'):strfind(sub_contents_tck(1).folder, 'sub')+6);
    
    % Read in all the .tck files available for this subject.
    for j = 1:length(sub_contents_tck)
        
        % Read in data for this subject and this tck.
        track_temp = read_mrtrix_tracks([sub_contents_tck(j).folder filesep sub_contents_tck(j).name]);
        
        % Merge.
        if j == 1
            
            mergedtck = track_temp.data;
            
        else
            
            mergedtck = [mergedtck track_temp.data];
            
        end
        
    end
    
    % Replace the streamline data in the track_temp struct so
    % that we keep the header information. (NOTE: This will
    % need to eventually be dealt with so that the header
    % information is accurate for the merged tck file and not
    % just the last tck file that was added to the merged tck
    % file.)
    track_temp.data = mergedtck;
    
    % Update basic header information.
    track_temp.count = length(track_temp.data);
    track_temp.total_count = length(track_temp.data);
    
    % Write out merged tck file.
    write_mrtrix_tracks(track_temp, [rootDir 'mergedtcks/' subID '_track.tck'])
    
    % Restart. 
    clear track_temp mergedtck sub_contents_tck subID
    
end


