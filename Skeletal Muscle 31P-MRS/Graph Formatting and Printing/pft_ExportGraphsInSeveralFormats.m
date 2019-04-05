function pft_ExportGraphsInSeveralFormats(HandleToFigure, Folder, FileNameStub)

% The figure should be created and deleted 'outside', before and after the function call
Suffix = { 'png', 'tif', 'emf', 'eps',  'pdf', 'svg' };
Format = { 'png', 'tif', 'emf', 'epsc', 'pdf', 'svg' };

NFORMATS = length(Format);

wb = waitbar(0, 'Exporting graphs');

for n = 1:NFORMATS
  PathName = fullfile(Folder, strcat(FileNameStub, '.', Suffix{n}));
  saveas(HandleToFigure, PathName, Format{n});
  waitbar(double(n)/double(NFORMATS), wb, sprintf('%1d of %1d files exported', n, NFORMATS));
end

waitbar(1, wb, sprintf('%1d of %1d graphs exported', NFORMATS, NFORMATS));
pause(0.5);
delete(wb);

end

