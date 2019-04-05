%% Clear the workspace as usual

clear all
close all
clc

fclose('all');

%% Prompt for a top-level folder with studies nested underneath

User = getenv('Username');
Base = fullfile('C:', 'Users', User, 'Desktop');
Root = uigetdir(Base, 'Choose a top-level folder with studies inside');

if ~ischar(Root)
  h = msgbox('No folder selected', 'Exit', 'modal');
  uiwait(h);
  delete(h);
  return;
end

%% Create a master XLSX summary file with 2 tabs

XlsxSummaryPath = fullfile(Root, 'Exercise 31-P MRS Summary.xlsx');

if (exist(XlsxSummaryPath, 'file') ~= 2)
  MonoHead = { 'Scanner ID', 'Subject ID', ...
               'Steady-State PCr Value / AU', 'Depletion Percentage', ...
               'Rate / s^-1', 'Time Constant / s' };

  BiHead = { 'Scanner ID', 'Subject ID', ...
             'Steady-State PCr Value / AU', 'Depletion Percentage', ...
             'Slow Fraction / PC', 'Slow Rate / s^-1', 'Slow Time Constant / s', ...
             'Fast Fraction / PC', 'Fast Rate / s^-1', 'Fast Time Constant / s' };
         
  DCHead = { 'Scanner ID', 'Subject ID', ...
             'Mono-Exponential kPCr / s^-1', 'At Dynamic Number', 'Time / s', ...
             'Bi-Exponential kPCr / s^-1', 'At Dynamic Number', 'Time / s', ...
             'Breakdown of PCr / PC' };
       
  DCMeanHead = { 'Scanner ID', 'Subject ID', ...
                 'Mean Mono-Exponential kPCr / s^-1', 'Mean Bi-Exponential kPCr / s^-1' };
      
  xlswrite(XlsxSummaryPath, MonoHead, 'Mono-Exponential Fit - Bout 1');
  xlswrite(XlsxSummaryPath, MonoHead, 'Mono-Exponential Fit - Bout 2');

  xlswrite(XlsxSummaryPath, BiHead, 'Bi-Exponential Fit - Bout 1');
  xlswrite(XlsxSummaryPath, BiHead, 'Bi-Exponential Fit - Bout 2');

  xlswrite(XlsxSummaryPath, DCHead, 'DC Statistics - Bout 1');
  xlswrite(XlsxSummaryPath, DCHead, 'DC Statistics - Bout 2');

  xlswrite(XlsxSummaryPath, DCMeanHead, 'DC - Mean Statistics');
end

%% List the sub-folders under the main working folder

Listing = dir(Root);

Entries = { Listing.name };
Folders = [ Listing.isdir ];
Entries = Entries(Folders);
Entries = Entries';
Entries = sort(Entries);

SingleDot = strcmpi(Entries, '.');
Entries(SingleDot) = [];
DoubleDot = strcmpi(Entries, '..');
Entries(DoubleDot) = [];

SubDirs = Entries;

NFOLDERS = numel(SubDirs);

%% Process the individual folders, allowing for the possibility of a deeply-nested failure with a try-catch block

[ Num, Txt, Raw ] = xlsread('Genscan II Subjects Directory.xlsx', 'Genscan II Subjects');

Head = Raw(1, :);
Data = Raw(2:end, :);

PrismaNumbers = Data(:, 1);
BRUNS         = Data(:, 2);

PrismaNumbers = cellfun(@strtrim, PrismaNumbers, 'UniformOutput', false);   % Trim possible leading and trailing spaces from the Excel entries
BRUNS         = cellfun(@strtrim, BRUNS, 'UniformOutput', false);           % Trim possible leading and trailing spaces from the Excel entries

for n = 1:NFOLDERS
  p = strfind(SubDirs{n}, '_');
  q = p(1) - 1;
  Leaf = SubDirs{n}(1:q);
  
  r = find(strcmpi(PrismaNumbers, Leaf), 1, 'first');
  
  ScannerID = PrismaNumbers{r};
  SubjectID = BRUNS{r};    
    
  [ MES1, BES1, MES2, BES2, DCStats1, DCStats2, DCMeanStats ] = pft_AnalyzeOneStudy(fullfile(Root, SubDirs{n}), ScannerID, SubjectID);  
    
  MonoData1 = horzcat({ ScannerID, SubjectID }, MES1);
  xlsappend(XlsxSummaryPath, MonoData1, 'Mono-Exponential Fit - Bout 1');
  MonoData2 = horzcat({ ScannerID, SubjectID }, MES2);
  xlsappend(XlsxSummaryPath, MonoData2, 'Mono-Exponential Fit - Bout 2');      
      
  BiData1 = horzcat({ ScannerID, SubjectID }, BES1);
  xlsappend(XlsxSummaryPath, BiData1, 'Bi-Exponential Fit - Bout 1');
  BiData2 = horzcat({ ScannerID, SubjectID }, BES2);
  xlsappend(XlsxSummaryPath, BiData2, 'Bi-Exponential Fit - Bout 2');    
    
  DCData1 = horzcat({ ScannerID, SubjectID }, DCStats1);
  xlsappend(XlsxSummaryPath, DCData1, 'DC Statistics - Bout 1');
  DCData2 = horzcat({ ScannerID, SubjectID }, DCStats2);
  xlsappend(XlsxSummaryPath, DCData2, 'DC Statistics - Bout 2');
   
  DCMeanData = horzcat({ ScannerID, SubjectID }, DCMeanStats);
  xlsappend(XlsxSummaryPath, DCMeanData, 'DC - Mean Statistics');
end

%% Signal completion

fprintf('All done !\n');