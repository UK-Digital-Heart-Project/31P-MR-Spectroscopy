% TITLE: Muscle 31P MR Spectroscopy Reconstruction Code SIEMENS
% PURPOSE: Uses Siemens data import scripts to import muscle 31P data
% from a DICOM file, performs post-processing, and saves the resulting 
% data in .txt format for jMRUI import.
%
% REQUIREMENTS: MATLAB Image Toolbox, MATLAB Signal Toolbox. Also requires 
% FID-A MATLAB scripts - 'SiemensCsaReadFid' and 'SiemensCsaParse' - to be 
% in the path (cite FID-A MRM paper), along with Chen's ACME scripts and 
% the Wiegers correlation method script.
% AUTHOR: Donnie Cameron
% DATE: 01/09/2017
% LAST UPDATED: 21/05/2018
%=============================================================================

tic
clearvars

%% SELECT AND IMPORT SIEMENS MRS FILES
% Select a single .IMA file from the series of 240 spectra and the script will
% automatically import the rest.
[ fileName, fileDir ] = uigetfile( '*.IMA' );
cd( fileDir );

ex_files = dir( [ '*' fileName( end - 51 : end - 45 ) '*.IMA' ] );
cell_fids = cell( numel( ex_files ), 1 );
for i = 1 : numel( ex_files )
    data = dicominfo( ex_files( i ).name );
    [ cell_fids{ i }, info ] = SiemensCsaReadFid( data, true );
end

fids = cell2mat( cell_fids );
fids = reshape( fids, [ size( cell_fids{ 1 }, 1 ), numel( ex_files ) ] );

%% LINE-BROADENING
lb = 10;                                         % Line-broadening factor
dw = info.csa.RealDwellTime * 1E-9;              % Dwell time from header
bw = 1 / dw;                                     % Bandwidth from dwell time 
ti = ( 0 : 1 : size( fids, 1 ) -1 ) .* dw;       % Time vector

gauss_ap = exp( - ( ti .* lb ) .^ 2 );        % Gaussian apodisation function

% Convolve data with Gaussian apodisation fn.
fids_ap = bsxfun( @times, fids, gauss_ap' );

%% ZERO-FILLING
% Use two-times zero-filling to improve spectral resolution
zero_fill = zeros( size( fids_ap, 1 ) .* 3, size( fids_ap, 2 ) ); 

fids_ap_zf = [ fids_ap; zero_fill ];

%% MANUAL PHASE CORRECTION
% Phase correction terms (zero-order phase will be automatically corrected)
zo_phase = 0;          % Zero-order phase correction in degrees.
begin_time = 0.00045;  % First-order phase correction in s.

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

% ACME automatic phase correction (cite Chen et al. 2002, JMR)
phc0arr = zeros( 1, sz( 2 ) );
phc1arr = zeros( 1, sz( 2 ) );
for i = 1 : sz( 2 )
    [ spec_ap_zf_phcorr( :, i ), phc0, phc1 ] = ...
        ACME( spec_ap_zf_phcorr( :, i ), [ 30, 1E-6 ] );
    % Using a tiny start value for the 1st-order phase effectively fixes it.
    phc0arr( i ) = phc0;
    phc1arr( i ) = phc1;
end

% Back to time domain for frequency alignment
fids_ap_zf_phcorr = ifft( spec_ap_zf_phcorr );

%% AUTOMATIC FREQUENCY ALIGNMENT
% Based on the correlation method of Wiegers et al. (MAGMA 2017). 
% Initialise arrays for storing corrected data.
fids_ap_zf_ph_fcorr = zeros( size( fids_ap_zf_phcorr ) );
fids_ap_zf_ph_fcorr( :, 1 ) = fids_ap_zf_phcorr( :, 1 );

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
        fids_ap_zf_phcorr( :, 1 ), ...
        fids_ap_zf_phcorr( :, i ), ...
        t_axis ), ...
        delt_sv );
    
    % Apply frequency and phase correction terms to FIDs.
    fids_ap_zf_ph_fcorr( :, i ) = ...
    fids_ap_zf_phcorr( :, i ) .* ...
    exp( 1i * 2 * pi * t_axis' * delt{ i - 1 }( 1 ) ) .* ...
    exp( 1i * 2 * pi * delt{ i - 1 }( 2 ) / 360 );

end

spec_ap_zf_ph_fcorr = fft( fids_ap_zf_ph_fcorr);

%% Check phase correction & frequency alignment
% count = 1;
% figure
% for i = 1 : 8 : size( spec_ap_zf_phcorr, 2 )
%     subplot( 5, 6, count ); plot( real( spec_ap_zf_phcorr( :, i ) ) );
%     xlim( [ 524, size( spec_ap_zf_phcorr, 1 ) - 325 ] )
%     ylim( [ min( imag( spec_ap_zf_phcorr( : ) ) * 0.3 ), ...
%         max( abs( spec_ap_zf_phcorr( : ) ) ) ] )
%     hold on
%     plot( real( spec_ap_zf_ph_fcorr( :, i ) ) );
%     count = count + 1;
%     hold off
% end

%% Filter spectral data in the indirect time domain
% First remove first 3 spectra - PCr signal is not at steady state
spec_ap_zf_ph_fcorr = spec_ap_zf_ph_fcorr( :, 4:end );
% Now, filter data.
spec_ap_zf_ph_fcorr_filt = itd_filter( spec_ap_zf_ph_fcorr, 2 );

%% Prepare data for export to jMRUI
% FT spectra prior to saving and split into real and imaginary components.
fids_real = real( fft( fftshift( spec_ap_zf_ph_fcorr_filt ) ) );
fids_imag = imag( fft( fftshift(spec_ap_zf_ph_fcorr_filt ) ) );

% Real and imaginary components of spectra.
spec_real = real( spec_ap_zf_ph_fcorr_filt );
spec_imag = imag( spec_ap_zf_ph_fcorr_filt );

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
        mr_ind + 14 : mr_ind + 29 ] ), '_31P_EX_recon.txt' ], 'wt' );
    date = fileName( mr_ind + 14 : mr_ind + 23 ); % Exam date.
else
    % Underscore in file name: set indices accordingly.
    fid = fopen( [ fileName( [ 1 : mr_ind mr_ind + 6 : mr_ind + 10 ...
        mr_ind + 16 : mr_ind + 31 ] ), '_31P_EX_recon.txt' ], 'wt' );
    date = fileName( mr_ind + 16 : mr_ind + 25 ); % Exam date.
end

% Then, write .txt file header.
fprintf( fid, '%s\n\n', 'jMRUI Data Textfile', ...
    'Filename: ex31P_recon.txt' );
fprintf( fid, '%s\n', ...
    [ 'PointsInDataset:', ' ', num2str( size( fids_real, 1 ) ) ], ...
    [ 'DatasetsInFile:', ' ', num2str( size( fids_real, 2 ) ) ], ...
    [ 'SamplingInterval: ', num2str( dw * 1000 ) ], ...
    'ZeroOrderPhase:0', ...
    'BeginTime: 0', ...
    [ 'TransmitterFrequency:', ' ', ...
    num2str( info.csa.ImagingFrequency ), 'E6' ], ...
    'MagneticField: 3E0', ...
    'TypeOfNucleus: 1E0', ...
    [ 'Name of Patient:', ' ', fileName( 1 : mr_ind ) ], ...
    [ 'Date of Experiment:', ' ', date ], ...
    'Spectrometer: Siemens_Prisma_3T', ...
    'AdditionalInfo: na' );
fprintf( fid, '%s', 'SignalNames: ' );	
for i = 1 : ( size( fids_real, 2 ) )
    fprintf( fid, '%s', [ 'P31_EX_', fileName( 1 : 6 ), '_', num2str( i ) ] );
    fprintf( fid, '%s', ';' );
end

% Write MRS data to text file.
fprintf( fid, '\n\n\n%s\n', 'Signal and FFT' );		
fprintf( fid, '%s\t', 'sig(real)', 'sig(imag)', 'fft(real)' );
fprintf( fid, '%s\n', 'fft(imag)' );
for i = 1 : size( fids_real, 2 )
    fprintf( fid, '%s', 'Signal number: ', num2str( i ), ' out of ', ...
        num2str( size( fids_real, 2 ) ), ' in file' );
    fprintf( fid, '\n' );
    for j = 1 : size( fids_real, 1 )
        fprintf( fid, '%.4E\t', fids_real( j, i ) );
        fprintf( fid, '%.4E\t', fids_imag( j, i ) );
        fprintf( fid, '%.4E\t', spec_real( j, i ) );
        fprintf( fid, '%.4E\n', spec_imag( j, i ) );
    end
end
fclose( fid );

fprintf( 'Processing complete!\n' )

toc
