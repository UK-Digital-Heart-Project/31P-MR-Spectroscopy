% TITLE: Cardiac 31P MRS Analysis Script
% PURPOSE: Reads in text files exported from jMRUI and arranges data into a
% format suitable for analysis. Corrected PCr:ATP ratios are then printed.
%
% REQUIREMENTS: MATLAB Optimization Toolbox.
% AUTHOR: Donnie Cameron
% DATE: 17/05/2016
% LAST UPDATED: 04/06/2018
%=============================================================================

tic
clearvars
close all

[ fileName, fileDir ] = uigetfile( 'Septal_31P_Spectrum.txt' );
cd( fileDir );

%% DEFINE CONSTANTS
% Sequence TR, and T1 of PCr and y-ATP - for relaxation correction.
tr = 13500;
t1pcr = 5800;  % From El-Sharkawy et al. (2009) MRM.
t1yatp = 3100; % From El-Sharkawy et al. (2009) MRM.

%% READ 31P DATA FROM TEXT FILE EXPORTED FROM jMRUI - AMARES
%===========================================================
P31_results = card_results_parse( fileName );

pcr_amp = P31_results.pcr_amp;
yatp1_amp = P31_results.yatp1_amp;
yatp2_amp = P31_results.yatp2_amp;
dpg2_amp = P31_results.dpg2_amp;
dpg3_amp = P31_results.dpg3_amp;

%% CALCULATE CORRECTED PCr:ATP AND CRSDs
%===========================================================
% Calculate blood correction.
blood_corr = 0.15 * ( P31_results.yatp1_amp + P31_results.yatp1_amp );

% Calculate blood-corrected PCr:ATP
bc_pcr_atp = P31_results.pcr_amp / ...
    ( P31_results.yatp1_amp + P31_results.yatp1_amp - blood_corr );

% Now apply relaxation correction (Handbook of MRS In Vivo)
fr = ( 1 - exp( - tr / t1yatp ) ) / ( 1 - exp( - tr / t1pcr ) ); 
rc_bc_pcr_atp = bc_pcr_atp * fr; 

% Calculate Cramer-Rao standard deviation (%) for each peak.
pcr_crsd = P31_results.pcr_amp_sd / P31_results.pcr_amp * 100;
yatp_crsd = P31_results.yatp1_amp_sd / P31_results.yatp1_amp * 100; 
dpg2_crsd = P31_results.dpg2_amp_sd / P31_results.dpg2_amp * 100; 
dpg3_crsd = P31_results.dpg3_amp_sd / P31_results.dpg3_amp * 100; 
% N.B. Each peak in a multiplets has the same CRSD. 

% Calculate PCr SNR.
snr_pcr = P31_results.pcr_amp / P31_results.noise;

% Now, print the results to a text file.
res = table( pcr_amp, ...
    pcr_crsd, ...
    yatp1_amp, ...
    yatp2_amp, ...
    yatp_crsd, ...
    dpg2_amp, ...
    dpg2_crsd, ...
    dpg3_amp, ...
    dpg3_crsd, ...
    bc_pcr_atp, ...
    rc_bc_pcr_atp, ...
    snr_pcr );

writetable( res, 'CARDIAC31P_RESULTS.csv' )

fprintf( 'Processing complete!\n' )

toc
