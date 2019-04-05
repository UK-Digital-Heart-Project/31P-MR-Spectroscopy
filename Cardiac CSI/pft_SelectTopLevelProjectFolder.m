%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Prompt for a top-level folder with subject studies underneath

UserName = getenv('UserName');

Start = fullfile('C:', 'Users', UserName, 'Desktop');
Title = 'Select a top-level project folder';

Root = uigetdir(Start, Title);

if ~ischar(Root)
  h = msgbox('No folder chosen', 'Quitting', 'modal');
  uiwait(h);
  delete(h);
  return;
end

%% Write the chosen path to a local text file

fid = fopen('Top-Level Project Folder.txt', 'wt');

fprintf(fid, '%s\n', Root);

fclose(fid);

%% Signal completion

h = msgbox('Top-level folder chosen', 'Success !', 'modal');
uiwait(h);
delete(h);


  