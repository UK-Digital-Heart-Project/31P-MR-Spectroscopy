% Clear the workspace

clear all
close all
clc

fclose('all');

% Nominate 2 base working directories

Source = 'C:\Users\ptokarcz\Desktop\31-P\Genscan II\Cardiac 31-P CSI';
Target = 'C:\Users\ptokarcz\Desktop\Stable CSI Workflow';

%% Copy entire data folders from the source to the target area - HCM

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'HCM');

Head = Raw(1, :);
Data = Raw(2:end, :);

ColA = find(strcmpi(Head, '3T Number'), 1, 'first');

SubjectID = Data(:, ColA);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Destination = fullfile(Target, 'HCM');

N = numel(SubjectID);

wb = waitbar(0, 'Copying HCM folders');

for n = 1:N
  [ Status, Message, MessageID ] = copyfile(fullfile(Source, SubjectID{n}), fullfile(Destination, SubjectID{n}));
  waitbar(double(n)/double(N), wb, sprintf('%1d of %1d folders copied', n, N));
end

waitbar(1, wb, sprintf('%1d of %1d folders copied', N, N'));

pause(1.0);

delete(wb);

%% Copy entire data folders from the source to the target area - Negative

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'Negative');

Head = Raw(1, :);
Data = Raw(2:end, :);

ColA = find(strcmpi(Head, '3T Number'), 1, 'first');

SubjectID = Data(:, ColA);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Destination = fullfile(Target, 'Negative');

N = numel(SubjectID);

wb = waitbar(0, 'Copying Negative folders');

for n = 1:N
  [ Status, Message, MessageID ] = copyfile(fullfile(Source, SubjectID{n}), fullfile(Destination, SubjectID{n}));
  waitbar(double(n)/double(N), wb, sprintf('%1d of %1d folders copied', n, N));
end

waitbar(1, wb, sprintf('%1d of %1d folders copied', N, N'));

pause(1.0);

delete(wb);

%% Copy entire data folders from the source to the target area - TTNtv (High PSI)

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'TTNtv (High PSI)');

Head = Raw(1, :);
Data = Raw(2:end, :);

ColA = find(strcmpi(Head, '3T Number'), 1, 'first');

SubjectID = Data(:, ColA);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Destination = fullfile(Target, 'TTNtv (High PSI)');

N = numel(SubjectID);

wb = waitbar(0, 'Copying TTNtv (High PSI) folders');

for n = 1:N
  [ Status, Message, MessageID ] = copyfile(fullfile(Source, SubjectID{n}), fullfile(Destination, SubjectID{n}));
  waitbar(double(n)/double(N), wb, sprintf('%1d of %1d folders copied', n, N));
end

waitbar(1, wb, sprintf('%1d of %1d folders copied', N, N'));

pause(1.0);

delete(wb);

%% Copy entire data folders from the source to the target area - TTNtv (PSI < 90)

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'TTNtv (PSI < 90)');

Head = Raw(1, :);
Data = Raw(2:end, :);

ColA = find(strcmpi(Head, '3T Number'), 1, 'first');

SubjectID = Data(:, ColA);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Destination = fullfile(Target, 'TTNtv (PSI_LT_90');

N = numel(SubjectID);

wb = waitbar(0, 'Copying TTNtv (PSI < 90) folders');

for n = 1:N
  [ Status, Message, MessageID ] = copyfile(fullfile(Source, SubjectID{n}), fullfile(Destination, SubjectID{n}));
  waitbar(double(n)/double(N), wb, sprintf('%1d of %1d folders copied', n, N));
end

waitbar(1, wb, sprintf('%1d of %1d folders copied', N, N'));

pause(1.0);

delete(wb);
