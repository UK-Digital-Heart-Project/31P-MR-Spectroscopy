function [ FirstFile, Container ] = pft_FindImaContainer(Folder)

% Set a default value for an early return - an empty array will yield true with isempty and false with ischar
FirstFile = [];
Container = [];

% List all the files and folders below the top level - don't attempt to use a wildcard at this stage
[ Files, Folders, Bytes ] = rdir(Folder);

% If no files of any kind are found, return with the default value
if isempty(Files)
  return;
end

% Count the files and search for the first with an IMA suffix; return the containing folder as soon as that is found - this overwrites the default value
NFILES = numel(Files);

for n = 1:NFILES
  [Path, Name, Extension ] = fileparts(Files{n});
  
  if strcmpi(Extension, '.IMA')
    FirstFile = Name;
    Container = Path;
    return;
  end
end

end
  
  





