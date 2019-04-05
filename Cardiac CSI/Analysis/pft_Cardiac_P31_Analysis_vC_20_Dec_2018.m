function pft_Cardiac_P31_Analysis_vC_20_Dec_2018(Home, Root, SubFolder, BRUN)
% TITLE: Cardiac 31P MRS Analysis Script
% PURPOSE: Reads in text files exported from jMRUI and arranges data into a
% format suitable for analysis. Corrected PCr:ATP ratios are then printed.
%
% REQUIREMENTS: MATLAB Optimization Toolbox.
% AUTHOR:       Donnie Cameron
% DATE:         17/05/2016
% LAST UPDATED: 04/06/2018
%
% RE-WRITTEN AS A FUNCTION TO REDUCE THE NEED FOR USER INTERACTION AND IMPROVE AUTOMATION
%
% AUTHOR:       Pawel Tokarczuk
% DATE:         15/11/2018
% LAST UPDATED: 20/12/2018
%
% Parameters:
%
% Home:         A base-level working directory, where a calling script is located
% Root:         A top-level folder with subject sub-folders (bearing 3 T numbers) underneath
% SubFolder:    A scanner ID, beginning with "3T"
% BRUN:         A patient ID, beginning with "14"

%% INITIALIZE THE WORKING DIRECTORY AND TEXT-MODE PRE-PROCESSED SEPTAL SPECTRUM DATA FILE

Away = fullfile(Root, SubFolder);

cd(Away);

SourceFileName = 'Septal_31P_Spectrum.txt';

if (exist(SourceFileName, 'file') ~= 2)
  hMsgBox = pft_MsgBox('No results file found', 'Exit', 'modal');
  uiwait(hMsgBox);
  delete(hMsgBox);
  cd(Home);
  return;
end

%% DEFINE CONSTANTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence TR, and T1 of PCr and y-ATP - for relaxation correction.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr     = 13500.0;
t1pcr  = 5800.0;      % From El-Sharkawy et al. (2009) MRM.
t1yatp = 3100.0;      % From El-Sharkawy et al. (2009) MRM.

%% READ 31P DATA FROM TEXT FILE EXPORTED FROM jMRUI - AMARES

fileDir  = Away;
fileName = 'Septal_31P_Spectrum.txt';

P31_results = card_results_parse( fullfile(fileDir, fileName) );

pcr_amp   = P31_results.pcr_amp;
yatp1_amp = P31_results.yatp1_amp;
yatp2_amp = P31_results.yatp2_amp;
dpg2_amp  = P31_results.dpg2_amp;
dpg3_amp  = P31_results.dpg3_amp;

%% CALCULATE CORRECTED PCr:ATP AND CRSD's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the raw PCr:ATP ratio                                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcr_atp = P31_results.pcr_amp / ( P31_results.yatp1_amp + P31_results.yatp2_amp );                  % PFT - 20/12/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the blood correction.                                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blood_corr = 0.15 * ( P31_results.dpg2_amp + P31_results.dpg3_amp );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate blood-corrected PCr:ATP.                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bc_pcr_atp = P31_results.pcr_amp / ( P31_results.yatp1_amp + P31_results.yatp1_amp - blood_corr ); % Duplicate terms (1, 2): typo ?
  bc_pcr_atp = P31_results.pcr_amp / ( P31_results.yatp1_amp + P31_results.yatp2_amp - blood_corr ); % PFT - 20/12/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Now apply relaxation correction (Handbook of MRS In Vivo).                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fr = ( 1 - exp( - tr / t1yatp ) ) / ( 1 - exp( - tr / t1pcr ) ); 
rc_bc_pcr_atp = bc_pcr_atp * fr; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Cramer-Rao standard deviation (%) for each peak.                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcr_crsd  = P31_results.pcr_amp_sd / P31_results.pcr_amp * 100;
yatp_crsd = P31_results.yatp1_amp_sd / P31_results.yatp1_amp * 100; 
dpg2_crsd = P31_results.dpg2_amp_sd / P31_results.dpg2_amp * 100; 
dpg3_crsd = P31_results.dpg3_amp_sd / P31_results.dpg3_amp * 100; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  N.B. Each peak in a multiplets has the same CRSD.                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate PCr SNR.                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr_pcr = P31_results.pcr_amp / P31_results.noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Now, print the results to a text file and a parallel XLSX sheet.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results = table( pcr_amp, ...
                 pcr_crsd, ...
                 yatp1_amp, ...
                 yatp2_amp, ...
                 yatp_crsd, ...
                 dpg2_amp, ...
                 dpg2_crsd, ...
                 dpg3_amp, ...
                 dpg3_crsd, ...
                 pcr_atp, ...
                 bc_pcr_atp, ...
                 rc_bc_pcr_atp, ...
                 snr_pcr );
         
OutputFileName = sprintf('%s - %s - CARDIAC_31P_RESULTS.csv', SubFolder, BRUN);

writetable( Results, OutputFileName );

Head = { '3T Number', 'BRUN', ...
         'PCr Amplitude', 'PCr C-R S.D. (per cent)', ...
         'Gamma-ATP-1 Amplitude', 'Gamma-ATP-2 Amplitude', 'Gamma-ATP C-R S.D. (per cent)', ...
         'DPG-2 Amplitude', 'DPG-2 C-R S.D. (per cent)', 'DPG-3 Amplitude', 'DPG-3 C-R S.D. (per cent)', ...
         'Raw PCr/ATP Ratio', 'Blood-Corrected PCr/ATP Ratio', 'Relaxation and Blood-Corrected PCr/ATP Ratio', ...
         'PCr SNR' };
     
Data = table2cell(Results);
Data = horzcat({ SubFolder, BRUN }, Data);

Full = vertcat(Head, Data);

OutputFileName = sprintf('%s - %s - CARDIAC_31P_RESULTS.xlsx', SubFolder, BRUN);

xlswrite(OutputFileName, Full, 'Cardiac CSI Results');

xlsappend(fullfile(Root, 'Collated Cardiac CSI Results.xlsx'), Data, 'Cardiac CSI Results');

%% SIGNAL COMPLETION

fprintf( 'Processing complete!\n' );

cd(Home);

end
