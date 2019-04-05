%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Read in all the patient genotype data

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'GenScan II Subjects');

Head = Raw(1, :);
Data = Raw(2:end, :);

%% Extract the extended genotype from the header and group the data accordingly

ColA = find(strcmpi(Head, 'Genotype'), 1, 'last');

ExtendedGenotypes = Data(:, ColA);

[ UniqueGenotypes, IA, IC ] = unique(ExtendedGenotypes);

NTYPES = numel(UniqueGenotypes);

for N = 1:NTYPES
  Part = Data(IC == N, :);
  Name = UniqueGenotypes{N};
  Full = vertcat(Head, Part);
  xlswrite('GenScan II Subjects Directory with Genotype - Extended.xlsx', Full, Name);
end

%% Signal completion

fprintf('All done !\n');




  