% TITLE: Muscle 31P MRS Post-Processing and Analysis Script
% PURPOSE: Reads in text files exported from jMRUI and arranges data into a
% format suitable for analysis. QC figures are generated, curve-fitting is
% performed, and kPCr and other measures are printed.
%
% REQUIREMENTS: MATLAB Optimization Toolbox.
% AUTHOR: Donnie Cameron
% DATE: 17/05/2016
% LAST UPDATED: 16/04/2018
%=============================================================================

tic
clearvars
close all

%****************************** USER INPUT ***********************************
% ENTER "t_dyn", IN SECONDS, HERE: i.e. NO. OF EXCITATIONS * REPETITION TIME
t_dyn = 3;
%*****************************************************************************

[ fileName, fileDir ] = uigetfile( '*RESULTS.txt' );
cd( fileDir );

%% READ 31P DATA FROM TEXT FILE EXPORTED FROM jMRUI - AMARES
%===========================================================
[ P31_results, n_dyn ] = results_parse( fileName );

%% SPLIT DATA INTO TWO REST-EXERCISE-REST BOUTS
%===========================================================
% Get number of rows in results table.
res_rows = size( P31_results, 1 );

% Split table in two bouts for separate automatic fitting.
P31_bout = cell( 1, 2 );
P31_bout{ 1 } = P31_results( 1 : 117, : ); % First three pts removed in recon.
P31_bout{ 2 } = P31_results( 118 : end, : );

%% PLOT TOTAL 31P TIME COURSE
%=============================
% Need data from individual peaks ... ADD THIS FEATURE
figure( 'Name', '31P Metabolite Time Course', 'NumberTitle', 'off' )
scatter( 1 : n_dyn, P31_results.pcr_amp, ...
    'r', 'filled', ...
    'MarkerFaceColor', 'r' );
hold on
scatter( 1 : n_dyn, P31_results.pi_amp, ...
    'b', 'filled', ...
    'MarkerFaceColor', 'b' );
xlabel( 'Time (s)' )
ylabel( 'Signal Intensity (a.u.)' )
ylim( [ 0 1.05 .* max( P31_results.pcr_amp ) ] )
xlim( [ 1 n_dyn ] )
legend( 'PCr', 'Pi' , 'Location', 'southeast' )
hold off

%% CURVE-FITTING PCr RECOVERY
%=============================
% Curve-fitting for each bout.
KPCr_results = cell( 1, numel( P31_bout ) );

fid = fopen( 'SUMMARY_RESULTS.txt', 'wt' ); 
for i = 1 : numel( P31_bout )
    
    % Determine minimum of PCr amplitude.
    [ pcr_min, pcr_min_index ] = min( P31_bout{ i }.pcr_amp );
    % Determine the 'rest' PCr amplitude and normalise data.
    pcr_rest = median( P31_bout{ i }( 1 : 15, : ).pcr_amp );
    % Initialise arrays for curve-fitting.
    pcr_init = P31_bout{ i }.pcr_amp( 1 : end ); % Normalise? / pcr_rest;
    t_arr_init = ones( numel( pcr_init ), 1 );
    t_arr_init = find( t_arr_init ) .* t_dyn;
    
    fprintf( fid, [ '\nBout ', num2str( i ), '\n' ] );
    
    [ KPCr_results{ i }, fits ] = p31ex_curvefit( P31_bout{ i }.pcr_amp, ...
        pcr_init, pcr_min_index, t_arr_init, t_dyn );  
    
    fprintf( fid, [ fits{ 1 } fits{ 2 } ] );
    
    fprintf( fid, [ 'Breakdown of PCr = ', ...
        num2str( ( pcr_rest - pcr_min ) / pcr_rest * 100 ), '%%\n' ] );
    
    fprintf( fid, [ 'Minimum pH = ', ...
        num2str( min( P31_bout{ i }.pH_dyn ) ), '\n' ] );
    
end

fprintf( fid, [ '\nMean mono kPCr = ', num2str( ...
    ( KPCr_results{ 1 }.mono_KPCr + ...
    KPCr_results{ 2 }.mono_KPCr ) / 2 ), '\n' ] );
fprintf( fid, [ 'Mean bi kPCr = ', num2str( ...
    ( KPCr_results{ 1 }.bi_KPCr + ...
    KPCr_results{ 2 }.bi_KPCr ) / 2 ), '\n\n' ] );
fclose( fid );

fprintf( 'Processing complete!\n' )

toc
