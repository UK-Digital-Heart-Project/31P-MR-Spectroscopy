%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Define some processing parameters

%  Here - a place to which to return (there are many migrations between folders in the worker functions)
Home = pwd;

%  A top-level project folder, which may already be designated in a local file - otherwise, it needs to be selected and noted down
if (exist('Top-Level Project Folder.txt', 'file') ~= 2)
  pft_SelectTopLevelProjectFolder;
end

if (exist('Top-Level Project Folder.txt', 'file') ~= 2)
  h = msgbox('No top-level folder selected', 'Quitting', 'modal');   
  uiwait(h);
  delete(h);
  return;
end

fid = fopen('Top-Level Project Folder.txt', 'rt');
    
Root = fgetl(fid);

fclose(fid);

%  A subject sub-folder
SubFolderFullPath = uigetdir(Root, 'Select a subject sub-folder');

if ~ischar(SubFolderFullPath)
  h = msgbox('No sub-folder selected', 'Quitting', 'modal');
  uiwait(h);
  delete(h);
  return;
end

p = strfind(SubFolderFullPath, filesep);
q = p(end) + 1;

SubFolder = SubFolderFullPath(q:end);

%  Finally, a BRUN for documentation

[ Num, Txt, Raw ] = xlsread('Genscan II Subjects Directory.xlsx', 'Genscan II Subjects');

Head = Raw(1, :);
Data = Raw(2:end, :);

ScannerID = Data(:, 1);
SubjectID = Data(:, 2);

ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

r = find(strcmpi(ScannerID, SubFolder), 1, 'first');

BRUN = SubjectID{r};

%% Now perform the CSI pre-processing step

pft_Cardiac_P31_Recon_Siemens_vC_20_Dec_2018(Home, Root, SubFolder, BRUN);



