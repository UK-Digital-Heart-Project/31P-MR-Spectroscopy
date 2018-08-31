% TITLE: Muscle 31P MRS Combined Recon. & Analysis Code SIEMENS
% PURPOSE: Uses Siemens data import scripts to import muscle 31P data
% from a DICOM file, performs post-processing, fits the resulting 
% phosphocreatine time courses, and saves fit parameters in a text file.
%
% REQUIREMENTS: MATLAB Image Toolbox, MATLAB Signal Toolbox. Also requires 
% FID-A MATLAB scripts - 'SiemensCsaReadFid' and 'SiemensCsaParse' - to be 
% in the path (cite FID-A MRM paper), along with the ACME automatic phase
% correction and Wiegers correlation method scripts.
%
% AUTHORS: Donnie Cameron & Pawel Tokarczuk
% DATE CREATED: 01/09/2017
% LAST UPDATED: 31/08/2018
%=============================================================================

tic
clearvars
close all

%****************************** USER INPUT ***********************************
% ENTER FIRST-ORDER PHASE CORRECTION TERM, "begin_time" HERE
begin_time = 0.00045;  % First-order phase correction in s.

% Windows for phosphocreatine and inorganic phosphate peak integration
pcr_win = 955 : 1015;
pi_win = 790 : 850;

% ENTER "t_dyn", IN SECONDS, HERE: i.e. NO. OF EXCITATIONS * REPETITION TIME
t_dyn = 3;
%*****************************************************************************

%% SELECT AND IMPORT SIEMENS MRS FILES
% Select a single .IMA file from the series of 240 spectra and the script will
% automatically import the rest.
[ fileName, fileDir ] = uigetfile( '*.IMA', 'Select a Siemens .IMA File' );
cd( fileDir );

disp( 'Importing Siemens 31P MRS data ...' );
ex_files = dir( [ '*' fileName( end - 51 : end - 45 ) '*.IMA' ] );
cell_fids = cell( numel( ex_files ), 1 );
for i = 1 : numel( ex_files )
    data = SiemensCsaParse( ex_files( i ).name );
    [ cell_fids{ i }, info ] = SiemensCsaReadFid( data, true );
end

fids = cell2mat( cell_fids );
fids = reshape( fids, [ size( cell_fids{ 1 }, 1 ), numel( ex_files ) ] );

%% DC CORRECTION
% Corrects for any DC offset in the time-domain data.
disp( 'DC-correcting exercise 31P MRS data ...' );

dc_off = mean( real( fids( end - 63 : end, : ) ) ); % Use last 64pts for offset
fids_dc = zeros( size( fids ) );
for i = 1 : size( fids, 2 )
    fids_dc( :, i ) = fids( :, i ) - dc_off( i );   % Subtract offset from FIDs
end

%% ZERO-FILLING
% Use two-times zero-filling to improve spectral resolution
disp( 'Zero-filling exercise 31P MRS data ...' );

zero_fill = zeros( size( fids, 1 ) .* 3, size( fids, 2 ) ); 
fids_zf = [ fids; zero_fill ];

%% MANUAL FIRST-ORDER PHASE CORRECTION
% First-order phase correction from user input, above.
disp( 'Performing first-order phase-correction from user input ...' );

% Calculate frequency axis.
dw = info.csa.RealDwellTime * 1E-9;            % Dwell time from header
bw = 1 / dw;                                   % Spectral bandwidth
sz = size( fids_zf );
f = ( - bw / 2 ) + ( bw /( 2 * sz( 1 ) ) ) : bw / ( sz( 1 ) ) : ...
    ( bw / 2 ) - ( bw / ( 2 * sz( 1 ) ) );
f_axis_ppm = f / info.csa.ImagingFrequency;
f_axis_ppm = f_axis_ppm + 1.2; % Position of PCr reference peak in ppm.
f_axis_hz = ( f_axis_ppm - 1.2 ) * info.csa.ImagingFrequency;

% Apply first-order phase correction to spectra
spec_zf_phcorr1 = bsxfun( @times, ...
    fftshift( fft( fids_zf ), 1 ), ...
    ones( size( fids_zf, 1 ), 1 ) .* ...
    exp( - 1i * 2 * pi * f_axis_hz' * begin_time ) );

%% AUTOMATIC ZERO-ORDER PHASE-CORRECTION
% Based on the ACME phase correction algorithm (cite Chen et al. 2002, JMR)
disp( 'Performing automatic zero-order phase-correction ...' );

spec_zf_phcorr2 = zeros( size( spec_zf_phcorr1 ) );

% Phase-correct first spectrum only, and apply correction to all.
[ spec_temp, phc0, phc1 ] = ...
    ACME( spec_zf_phcorr1( pcr_win, 1 ), [ 0, 1E-20 ] );
% Using a tiny start value for the 1st-order phase effectively fixes it.

% Apply phase correction to all spectra.
for i = 1 : sz( 2 )      
    spec_zf_phcorr2( :, i ) = bsxfun( @times, ...
        spec_zf_phcorr1( :, i ), ...
        ( ones( size( spec_zf_phcorr2, 1 ), 1 ) * ... 
        exp( - 1i * 2 * pi * phc0 / 360 ) ) );
end

% Back to time domain for frequency alignment
fids_zf_phcorr = ifft( spec_zf_phcorr2 );

%% AUTOMATIC FREQUENCY ALIGNMENT
% Based on the correlation method of Wiegers et al. (MAGMA 2017). 
disp( 'Performing automatic frequency alignment ...' );

% Initialise arrays for storing corrected data.
fids_zf_ph_fcorr = zeros( size( fids_zf_phcorr ) );
fids_zf_ph_fcorr( :, 1 ) = fids_zf_phcorr( :, 1 );

% Calculate time axis.
t_axis = 0 : dw : ( sz( 1 ) - 1 ) * dw;

delt = cell( sz( 2 ) - 1, 1 ); % Save calculated terms in cell array.
delt_sv = [ 0, 0 ]; % Starting values for delta freq. and delta phase.

% Use the first spectrum in the series as a reference. Test correlation of
% all subsequent spectra versus the first.
for i = 2 : sz( 2 )
    
    % Determine frequency and phase shifts by maximising the cor_fun
    % function (fminsearch where function has negative sign).
    delt{ i - 1 } = fminsearch( @( delt_sv )cor_fun( delt_sv, ...
        fids_zf_phcorr( :, 1 ), ...
        fids_zf_phcorr( :, i ), ...
        t_axis ), ...
        delt_sv );
    
    % Update starting values based on previous result.  
    delt_sv = delt{ i - 1 };
    
    % Apply frequency and phase correction terms to FIDs.
    fids_zf_ph_fcorr( :, i ) = ...
    fids_zf_phcorr( :, i ) .* ...
    exp( 1i * 2 * pi * t_axis' * delt{ i - 1 }( 1 ) ) .* ...
    exp( 1i * 2 * pi * delt{ i - 1 }( 2 ) / 360 );

end

spec_zf_ph_fcorr = fft( fids_zf_ph_fcorr );

%% INDIRECT TIME DOMAIN FILTER
% Filters spectral data in the exercise domain
disp( 'Filtering exercise 31P MRS data in the indirect time domain ...' );

% First remove first 3 spectra - PCr signal is not at steady state
spec_zf_ph_fcorr = spec_zf_ph_fcorr( :, 4 : end );
% Now, filter data.
spec_zf_ph_fcorr_filt = itd_filter( spec_zf_ph_fcorr, 1.8 );

%% PEAK INTEGRATION
% Integrates under PCr and Pi peaks for each spectrum
disp( 'Integrating peak areas ...' );

% Real components of spectra.
spec_real = real( spec_zf_ph_fcorr_filt );

% Integrate under PCr and Pi peaks
PCr_amp = zeros( 1, size( spec_real, 2 ) );
Pi_amp = zeros( 1, size( spec_real, 2 ) );
for i = 1 : size( spec_real, 2 )
    PCr_amp( i ) = trapz( spec_real( pcr_win, i ) );
    Pi_amp( i ) = trapz( spec_real( pi_win, i ) );
end

%% SPLIT DATA INTO TWO REST-EXERCISE-REST BOUTS
% Splits data into two bouts for separate automatic fitting.
disp( 'Splitting data into two bouts ...' );

% Get number of rows in results table.
n_dyn = size( PCr_amp, 2 );

% Split table in two for separate fitting.
P31_bout = cell( 1, 2 );
P31_bout{ 1 }.pcr_amp = PCr_amp( 1, 1 : 117 )'; % First 3 pts removed in recon.
P31_bout{ 1 }.pi_amp = Pi_amp( 1, 1 : 117 )';
P31_bout{ 2 }.pcr_amp = PCr_amp( 1, 118 : end )';
P31_bout{ 2 }.pi_amp = Pi_amp( 1, 118 : end )';

%% PLOT TOTAL 31P TIME COURSE
% Plots the complete time course of the 31P MRS data - both exercise bouts.
disp( 'Plotting complete time course of exercise 31P MRS data ...' );

% Create new figure
figure( 'Name', '31P Metabolite Time Course', 'NumberTitle', 'off' )
scatter( 1 : n_dyn, PCr_amp, ...
    'r', 'filled', ...
    'MarkerFaceColor', 'r' );
hold on
scatter( 1 : n_dyn, Pi_amp, ...
    'b', 'filled', ...
    'MarkerFaceColor', 'b' );
xlabel( 'Time (s)' )
ylabel( 'Signal Intensity (a.u.)' )
ylim( [ 0 1.05 .* max( PCr_amp ) ] )
xlim( [ 1 n_dyn ] )
legend( 'PCr', 'Pi', 'Location', 'east' )
hold off

%% CURVE-FITTING PCr RECOVERY
% Fits PCr recovery after exercise with mono- and biexponential functions.
disp( 'Fitting PCr recovery curves ...' );

% Curve-fitting for each bout.
KPCr_results = cell( 1, numel( P31_bout ) );

fid = fopen( 'SUMMARY_RESULTS.txt', 'wt' ); 
for i = 1 : numel( P31_bout )
    
    bout_str = [ '... for bout ', num2str( i ) ];
    disp( bout_str )
    
    % Determine minimum of PCr amplitude.
    [ pcr_min, pcr_min_index ] = min( P31_bout{ i }.pcr_amp );
    % Determine the 'rest' PCr amplitude and normalise data.
    pcr_rest = median( P31_bout{ i }.pcr_amp( 1 : 15 ) );
    % Initialise arrays for curve-fitting.
    pcr_init = P31_bout{ i }.pcr_amp( 1 : end ); % Normalise? / pcr_rest;
    t_arr_init = ones( numel( pcr_init ), 1 );
    t_arr_init = find( t_arr_init ) .* t_dyn;
    
    fprintf( fid, [ '\nBout ', num2str( i ), '\n' ] );
    
    [ KPCr_results{ i }, fits ] = p31ex_curvefit( P31_bout{ i }.pcr_amp, ...
        pcr_init, pcr_min_index, t_arr_init, t_dyn, i );  
    
    fprintf( fid, [ fits{ 1 } fits{ 2 } ] );
    
    fprintf( fid, [ 'Breakdown of PCr = ', ...
        num2str( ( pcr_rest - pcr_min ) / pcr_rest * 100 ), '%%\n' ] );
    
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
