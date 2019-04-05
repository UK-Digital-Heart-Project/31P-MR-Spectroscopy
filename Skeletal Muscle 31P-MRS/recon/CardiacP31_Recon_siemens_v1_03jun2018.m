% TITLE: Cardiac 31P MR Spectroscopy Reconstruction Code SIEMENS
% PURPOSE: Uses Siemens data import scripts to import cardiac 31P data
% from a DICOM file, performs post-processing, and saves the resulting 
% data in .txt format for jMRUI import. 
%
% REQUIREMENTS: MATLAB Image Toolbox, MATLAB Signal Toolbox. Also requires 
% FID-A MATLAB scripts - 'SiemensCsaReadFid' and 'SiemensCsaParse' - to be 
% in the path (cite FID-A MRM paper), along with Chen's ACME scripts and 
% the Wiegers correlation method script.
%
% AUTHOR: Donnie Cameron
% DATE: 12/03/2018
% LAST UPDATED: 03/06/2018
%===============================================================================

tic
clearvars
close all

%% SELECT AND IMPORT SIEMENS MRS FILES
% Select a 2D-CSI .IMA file
[ fileName, fileDir ] = uigetfile( '*.IMA' );
cd( fileDir);

data = dicominfo( fileName );
[ fids, info ] = SiemensCsaReadFid( data, true );

% FIDs are in [ pts, cols, rows ] order, which is consistent with SpectrIm
% and the Siemens spectroscopy interface. 

%% PROMPT USER FOR SEPTAL VOXEL COORDINATES
% Provide a field for the user to enter 2D voxel coordinates.
prompt = { [ 'Enter voxel coordinates in [ x, y ] format. Select multiple' ...
    ' voxels by entering more bracketed values: e.g. [ x, y ] [ x, y ]' ] };
ptitle = 'Septal Voxel Coordinates';
dims = [ 1 60 ];
usr_input = inputdlg( prompt, ptitle, dims );
b1 = regexp( usr_input, '[' );
b2 = regexp( usr_input, ']' );

% Throw an error if there are no square brackets in the entered string, or
% if brackets are not paired.
if numel( b1{ 1 } ) == 0 && numel( b2{ 1 } ) == 0
    error( 'Please enclose coordinates in square brackets - [ ]' )
elseif numel( b1{ 1 } ) ~= numel( b2{ 1 } )
    error( 'Missing opening/closing bracket in entered coordinates' )
elseif any( b2{ 1 } - b1{ 1 } < 0 )
    error( 'Missing opening/closing bracket in entered coordinates' )
end

% Extract voxel coordinates from entered string. 
coord = cell( 1, numel( b1{ 1 } ) ); 
for i = 1 : numel( b1{ 1 } )
    coord{ i } = str2num( usr_input{ 1 }( b1{ 1 }( i ) : b2{ 1 }( i ) ) );
end

%% PRE-PROCESS SPECTRA
%=====================
% Pre-processing consists of decimation, CSI reconstruction, line-broadening, 
% zero-filling and phase-correction.

%% DECIMATE OVER-SAMPLED DATA
% Decimate data to prescribed number of sample points by averaging adjacent
% vector elements.
dec_factor = info.csa.DataPointColumns / ...
    info.csa.SpectroscopyAcquisitionDataColumns; % Acq. pts / user def. pts
if dec_factor == 2
    fids_dec = zeros( size( fids, 1 ) / dec_factor, ...
        size( fids, 2 ), ...
        size( fids, 3 ) );
    for i = 1 : size( fids, 2 )
        for j = 1 : size( fids, 3 )
            fids_dec( :, i, j ) = mean( ...
                [ fids( 1 : 2 : end, i, j ) fids( 2 : 2 : end, i, j ) ], 2 );                                
        end
    end
elseif dec_factor < 1
    warning( 'Data inconsistent with DICOM header info.' )
else
    fids_dec = fids;
end

% Note, no 1/2-voxel spatial shift in this script.

%% 2DCSI RECONSTRUCTION
% Fourier transformation in two dimensions for spatial recon. 
fids_recon = fftshift( fft( fids_dec, [], 2 ), 2 );
fids_recon = fftshift( fft( fids_recon, [], 3 ), 3 );

%% LINE-BROADENING
lb = 20;                                         % Line-broadening factor
dw = info.csa.RealDwellTime * dec_factor * 1E-9; % Dwell time from header
bw = 1 / dw;                                     % Bandwidth from dwell time 
ti = ( 0 : 1 : size( fids_dec, 1 ) -1 ) .* dw;   % Time vector

gauss_ap = exp( - ( ti .* lb ) .^ 2 );        % Gaussian apodisation function
fids_ap = zeros( size( fids_recon ) );        % Initialise apodised data array

% Convolve data with Gaussian apodisation fn.
for i = 1 : size( fids_recon, 2 )
    for j = 1 : size( fids_recon, 3 )
        fids_ap( :, i, j ) = fids_recon( :, i, j ) .* gauss_ap';
    end
end

%% ZERO-FILLING
% Use two-times zero-filling to improve spectral resolution
zero_fill = zeros( 1, size( fids_ap, 1 ) .* 3 ); 
% 2x Zero-filling vector, to be concatenated with data
fids_ap_zf = zeros( size( zero_fill, 2 ) + size( fids_ap, 1 ), ...
    size( fids_ap, 2 ), size( fids_ap, 2 ) );

for i = 1 : size( fids_ap, 2 )
    for j = 1 : size( fids_ap, 3 )
        fids_ap_zf( :, i, j ) = [ fids_ap( :, i, j ); zero_fill' ];
    end
end

%% MANUAL FIRST-ORDER PHASE CORRECTION
% Phase correction terms (zero-order phase will be automatically corrected)
zo_phase = 0;          % Zero-order phase correction in degrees.
begin_time = 0.0026;     % First-order phase correction in s.

% Apply manual zero-order phase correction (if needed). 
fids_ap_zf_phcorr = bsxfun( @times, fids_ap_zf, ...
    ( ones( size( fids_ap_zf, 1 ), 1 ) * ...
    exp( 1i * 2 * pi * zo_phase / 360 ) ) );

% Calculate frequency axis.
sz = size( fids_ap_zf );
f = ( - bw / 2 ) + ( bw /( 2 * sz( 1 ) ) ) : bw / ( sz( 1 ) ) : ...
    ( bw / 2 ) - ( bw / ( 2 * sz( 1 ) ) );
f_axis_ppm = f / info.csa.ImagingFrequency;
f_axis_ppm = f_axis_ppm + 1.2; % Position of PCr reference peak in ppm.
f_axis_hz = ( f_axis_ppm - 1.2 ) * info.csa.ImagingFrequency;

% Apply first-order phase correction to spectra
spec_ap_zf_phcorr = bsxfun( @times, ...
    fftshift( fft( fids_ap_zf_phcorr ), 1 ), ...
    ones( size( fids_ap_zf, 1 ), 1 ) .* ...
    exp( - 1i * 2 * pi * f_axis_hz' * begin_time ) );

%% ACME automatic phase correction (cite Chen et al. 2002, JMR)
phc0arr = zeros( sz( 2 ), sz( 3 ) );
phc1arr = zeros( sz( 2 ), sz( 3 ) );
for i = 1 : sz( 2 )
    for j = 1 : sz( 3 )
        [ spec_ap_zf_phcorr( :, i, j ), phc0, phc1 ] = ...
            ACME( spec_ap_zf_phcorr( :, i, j ), [ 30, 1E-6 ] );
        % Using tiny start val. for 1st-order phase effectively fixes it.
        phc0arr( i, j ) = phc0;
        phc1arr( i, j ) = phc1;
    end
end

% Back to time domain for frequency alignment
fids_ap_zf_phcorr = ifft( ifftshift( spec_ap_zf_phcorr, 1 ) );

%% Plot all 2D-CSI data for inspection
% count = 1;
% figure
% for y = 1:8
%     for x = 1:8 
%         subplot( 8, 8, count ); plot( real( spec_ap_zf_phcorr( :, x, y ) ) )
%         xlim( [ 524, size( spec_ap_zf_phcorr, 1 ) - 325 ] )
%         ylim( [ min( imag( spec_ap_zf_phcorr( : ) ) * 0.4 ), ...
%             max( abs( spec_ap_zf_phcorr( : ) ) ) ] )
%         count = count + 1;
%     end
% end

%% IDENTIFY SEPTAL F.I.D.s
% Take coordinates entered earlier and use them to select FIDs from 
% septal voxels.
fids_septum = zeros( size( fids_ap_zf_phcorr, 1 ), numel( coord ) );
% Loop over cell array of coords and create array of septal FIDs.
for i = 1 : numel( coord )
    fids_septum( :, i ) = ...
        fids_ap_zf_phcorr( :, coord{ i }( 1 ), coord{ i }( 2 ) );
end

%% AUTOMATIC FREQUENCY ALIGNMENT
% Based on the correlation method of Wiegers et al. (MAGMA 2017). 
% Initialise arrays for storing corrected data.
fids_septum_fcorr = zeros( size( fids_septum ) );

% Calculate time axis.
t_axis = 0 : dw : ( sz( 1 ) - 1 ) * dw;

% Save calculated terms in cell array.
delt = cell( size( fids_septum, 2 ), size( fids_septum, 3 ) ); 
delt_sv = [ 0, 0 ]; % Starting values for delta freq. and delta phase.

% Use the highest SNR spectrum in the series as a reference. Test correlation
% of all other spectra versus the reference.
for i = 1 : size( fids_septum, 2 )
    for j = 1 : size( fids_septum, 3 )
    
        [ m, I ] = max( max( abs( fids_septum ) ) ); % I = highest SNR fid
        
        % Determine frequency and phase shifts by maximising the cor_fun
        % function (fminsearch where function has negative sign).
        delt{ i, j } = fminsearch( @( delt_sv )cor_fun( delt_sv, ...
            fids_septum( :, I ), ...
            fids_septum( :, i, j ), ...
            t_axis ), ...
            delt_sv );

        % Apply frequency and phase correction terms to FIDs.
        fids_septum_fcorr( :, i, j ) = ...
        fids_septum( :, i, j ) .* ...
        exp( 1i * 2 * pi * t_axis' * delt{ i, j }( 1 ) ) .* ...
        exp( 1i * 2 * pi * delt{ i, j }( 2 ) / 360 );

    end
end

spec_septum_fcorr = fftshift( fft( fids_septum_fcorr), 1 );

%% AVERAGE SEPTAL SPECTRA
spec_septum_avg = mean( spec_septum_fcorr, 2 ); 

%% Plot selected septal voxels, and their average, for inspection
figure

subplot( 1, 2, 1 ); plot( real( spec_septum_fcorr ) )
xlim( [ 524, size( spec_septum_fcorr, 1 ) - 325 ] )
ylim( [ min( real( spec_septum_fcorr( : ) ) ), ...
    max( abs( spec_septum_fcorr( : ) ) ) ] )
title('Septal Spectra')
xlabel('Spectral Points') % x-axis label
ylabel('Signal Intensity (a.u)') % y-axis label

subplot( 1, 2, 2 ); plot( real( spec_septum_avg ) )
xlim( [ 524, size( spec_septum_avg, 1 ) - 325 ] )
ylim( [ min( real( spec_septum_avg( : ) ) ), ...
    max( abs( spec_septum_avg( : ) ) ) ] )
title('Average of Septal Spectra')
xlabel('Spectral Points') % x-axis label
ylabel('Signal Intensity (a.u)') % y-axis label

%% Prepare data for export to jMRUI
% FT spectra prior to saving and split into real and imaginary components.
fids_real = real( ifft( ifftshift( spec_septum_avg ) ) );
fids_imag = - imag( ifft( ifftshift( spec_septum_avg ) ) );
        
% Real and imaginary components of spectra.
spec_real = real( spec_septum_avg );
spec_imag = - imag( spec_septum_avg );

%% WRITE DATA TO .TXT FILE FOR jMRUI IMPORT
% First, generate .txt. file name using participant ID, series #, and date
% (in YYYY.MM.DD.HH.mm format) from DICOM file name.
% There are two naming conventions for .IMA files, depending on how files are
% exported. Test for these to ensure our filename comes out right.
mr_ind = regexp( fileName, '.MR._.' );            % Test for '.MR._.' string. 
if isempty( mr_ind )
    % No underscore in file name: set indices accordingly.
    mr_ind = regexp( fileName, '.MR.' );          % '.MR.' always in filename. 
    fid = fopen( [ fileName( [ 1 : mr_ind mr_ind + 4 : mr_ind + 8 ...
        mr_ind + 14 : mr_ind + 29 ] ), '_31P_card_recon.txt' ], 'wt' );
    date = fileName( mr_ind + 14 : mr_ind + 23 ); % Exam date.
else
    % Underscore in file name: set indices accordingly.
    fid = fopen( [ fileName( [ 1 : mr_ind mr_ind + 6 : mr_ind + 10 ...
        mr_ind + 16 : mr_ind + 31 ] ), '_31P_card_recon.txt' ], 'wt' );
    date = fileName( mr_ind + 16 : mr_ind + 25 ); % Exam date.
end

% Then, write .txt file header.
fprintf( fid, '%s\n\n', 'jMRUI Data Textfile', ...
    'Filename: card31P_recon.txt' );
fprintf( fid, '%s\n', ...
    [ 'PointsInDataset:', ' ', num2str( size( fids_real, 1 ) ) ], ...
    'DatasetsInFile: 1', ...
    [ 'SamplingInterval: ', num2str( dw * 1000 ) ], ...
    'ZeroOrderPhase: 0', ...
    'BeginTime: 0', ...
    [ 'TransmitterFrequency:', ' ', ...
    num2str( info.csa.ImagingFrequency ), 'E6' ], ...
    'MagneticField: 3E0', ...
    'TypeOfNucleus: 1E0', ...
    [ 'Name of Patient:', ' ', fileName( 1 : mr_ind ) ], ...
    [ 'Date of Experiment:', ' ', date ], ...
    'Spectrometer: Siemens_Prisma_3T', ...
    'AdditionalInfo: na' );
fprintf( fid, '%s', 'SignalNames: Septal_31P_Spectrum' );	

% Write MRS data to text file.
fprintf( fid, '\n\n\n%s\n', 'Signal and FFT' );		
fprintf( fid, '%s\t', 'sig(real)', 'sig(imag)', 'fft(real)' );
fprintf( fid, '%s\n', 'fft(imag)' );
fprintf( fid, '%s', 'Signal number: 1 out of 1 in file' );
fprintf( fid, '\n' );
for k = 1 : size( fids_real, 1 )
    fprintf( fid, '%.4E\t', fids_real( k, 1 ) );
    fprintf( fid, '%.4E\t', fids_imag( k, 1 ) );
    fprintf( fid, '%.4E\t', spec_real( k, 1 ) );
    fprintf( fid, '%.4E\n', spec_imag( k, 1 ) );
end
fclose( fid );

sprintf( 'Processing complete!\n' )

toc
