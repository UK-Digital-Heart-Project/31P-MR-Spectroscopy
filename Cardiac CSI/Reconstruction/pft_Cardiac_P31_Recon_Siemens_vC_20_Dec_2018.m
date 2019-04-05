function pft_Cardiac_P31_Recon_Siemens_vC_20_Dec_2018(Home, Root, SubFolder, BRUN)
% TITLE:        Cardiac 31P MR Spectroscopy Reconstruction Code SIEMENS
% PURPOSE:      Uses Siemens data import scripts to import cardiac 31P data from a DICOM file, 
%               performs post-processing, and saves the resulting data in .txt format for jMRUI import. 
%
% REQUIREMENTS: 
%
% MATLAB Image Toolbox, MATLAB Signal Toolbox. 
% Also requires FID-A MATLAB scripts - 'SiemensCsaReadFid' and 'SiemensCsaParse' - to be 
% in the path (cite FID-A MRM paper), along with Chen's ACME scripts and the Wiegers correlation method script.
%
% AUTHOR:       Donnie Cameron
% DATE:         12/03/2018
% LAST UPDATED: 03/06/2018
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

%% INITIALIZE THE WORKING DIRECTORY AND CSI FILE

Away = fullfile(Root, SubFolder);

cd(Away);

[ fileName, fileDir, filterIndex ] = uigetfile({'*.*';'*.dcm';'*.IMA'}, 'Select the CSI file (DCM or IMA)', Away);

if (filterIndex == 0)
  hMsgBox = pft_MsgBox('No file selected', 'Exit', 'modal');
  uiwait(hMsgBox);
  delete(hMsgBox);
  cd(Home);
  return;
end

%% READ IN THE CSI DATA FILE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FID's are in [ pts, cols, rows ] order, which is consistent with SpectrIm and the Siemens spectroscopy interface.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = dicominfo( fullfile(fileDir, fileName) );
[ fids, info ] = SiemensCsaReadFid( data, true );

%% PROMPT USER FOR SEPTAL VOXEL COORDINATES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Provide a field for the user to enter 2D voxel coordinates.                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prompt = sprintf([ 'Enter voxel coordinates in [ x, y ] format.\n', ... 
                   'Select multiple voxels by entering more bracketed values:\n', ... 
                   'e.g. [ x1, y1 ] [ x2, y2 ]' ]);

Title = 'Septal Voxel Coordinates';

Dims = [ 1 60 ];

usr_input = pft_InputDlg( Prompt, Title, Dims );
b1 = regexp( usr_input, '[' );
b2 = regexp( usr_input, ']' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Throw an error if there are no square brackets in the entered string, or if brackets are not paired.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel( b1{ 1 } ) == 0 && numel( b2{ 1 } ) == 0
    error( 'Please enclose coordinates in square brackets - [ ]' )
elseif numel( b1{ 1 } ) ~= numel( b2{ 1 } )
    error( 'Missing opening/closing bracket in entered coordinates' )
elseif any( b2{ 1 } - b1{ 1 } < 0 )
    error( 'Missing opening/closing bracket in entered coordinates' )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extract voxel coordinates from entered string.                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coord = cell( 1, numel( b1{ 1 } ) ); 
for i = 1 : numel( b1{ 1 } )
    coord{ i } = str2num( usr_input{ 1 }( b1{ 1 }( i ) : b2{ 1 }( i ) ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Write out the pixel co-ordinates for later use.                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputFileName = sprintf('%s - %s - Pixel Coordinates.txt', SubFolder, BRUN);

fid = fopen(OutputFileName, 'wt');

fprintf(fid, 'x/Col\ty/Row\n');

coord = coord(:);

NPTS = size(coord, 1);

for n = 1:NPTS
  fprintf(fid, '%d\t%d\n', coord{n}(1), coord{n}(2));
end

fclose(fid);

%% PRE-PROCESS SPECTRA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pre-processing consists of decimation, CSI reconstruction, line-broadening, zero-filling and phase-correction.                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize a waitbar to keep track of processing

wb = waitbar(0, sprintf('Pre-processing: %s - %s ...', SubFolder, BRUN));

NSTEPS = 17;
PHASES = 0;

%% DECIMATE OVER-SAMPLED DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Decimate data to prescribed number of sample points by averaging adjacent vector elements.                                       %
%  Note, there is no 1/2-voxel spatial shift in this script.                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dec_factor = info.csa.DataPointColumns / info.csa.SpectroscopyAcquisitionDataColumns; % Acquired points / user defined points
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

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Data decimated ...');

%% FID conjugation - to allow the use of an FT with the right sign and scaling for downstream analysis

fids_dec = conj(fids_dec);  % Donnie - is this correct ? Yes: 20/12/2018.
                            % Perhaps a matter of sign conventions between jMRUI/RDA/IMA/DCM ?

fids_recon = fids_dec;      % This is strictly redundant - one of the variables is left as a fossil,
                            % now that the spatial FT has been eliminated (following tests with a localization phantom)

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Data conjugated ...');

%% LINE-BROADENING

lb = 20;                                            % Line-broadening factor
dw = info.csa.RealDwellTime * dec_factor * 1E-9;    % Dwell time from header
bw = 1 / dw;                                        % Bandwidth from dwell time 
ti = ( 0 : 1 : size( fids_dec, 1 ) -1 ) .* dw;      % Time vector

gauss_ap = exp( - ( ti .* lb ) .^ 2 );              % Gaussian apodisation function
fids_ap = zeros( size( fids_recon ) );              % Initialise apodised data array

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Line-broadening initialized ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Convolve data with a Gaussian apodisation function.                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : size( fids_recon, 2 )
    for j = 1 : size( fids_recon, 3 )
        fids_ap( :, i, j ) = fids_recon( :, i, j ) .* gauss_ap';
    end
end

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Gaussian broadening applied ...');

%% ZERO-FILLING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Use two-times zero-filling to improve spectral resolution.                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zero_fill = zeros( 1, size( fids_ap, 1 ) .* 3 ); 

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Zero-filling - phase (1) completed ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2x zero-filling vector, to be concatenated with data.                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_ap_zf = zeros( size( zero_fill, 2 ) + size( fids_ap, 1 ), size( fids_ap, 2 ), size( fids_ap, 2 ) );

for i = 1 : size( fids_ap, 2 )
    for j = 1 : size( fids_ap, 3 )
        fids_ap_zf( :, i, j ) = [ fids_ap( :, i, j ); zero_fill' ];
    end
end

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Zero-filling - phase (2) completed ...');

%% MANUAL FIRST-ORDER PHASE CORRECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Phase correction terms (zero-order phase will be automatically corrected).                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zo_phase   = 0;                         % Zero-order phase correction in degrees.
begin_time = 0.0026;                    % First-order phase correction in seconds.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Apply manual zero-order phase correction (if needed).                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_ap_zf_phcorr = bsxfun( @times, fids_ap_zf, ( ones( size( fids_ap_zf, 1 ), 1 ) * exp( 1i * 2 * pi * zo_phase / 360 ) ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the frequency axis.                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size( fids_ap_zf );
np = sz(1);
f  = (-np/2:np/2-1)*bw/double(np);

f_axis_ppm = f / info.csa.ImagingFrequency;     % PCr appears at 1.2 ppm, both here and in jMRUI (before zero-referencing)
f_axis_hz  = f;                                 % Needed for the phase-correction steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Apply first-order phase correction to spectra.                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spec_ap_zf_phcorr = bsxfun( @times, ...
                            fftshift( fft( fids_ap_zf_phcorr ), 1 ), ...
                            ones( size( fids_ap_zf, 1 ), 1 ) .* exp( - 1i * 2 * pi * f_axis_hz' * begin_time ) );
                        
PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'First-order phase correction applied ...');

%% ACME automatic phase correction (cite Chen et al. 2002, JMR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Using a tiny start value for the 1st-order phase effectively fixes it.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phc0arr = zeros( sz( 2 ), sz( 3 ) );
phc1arr = zeros( sz( 2 ), sz( 3 ) );

for i = 1 : sz( 2 )
    for j = 1 : sz( 3 )
        [ spec_ap_zf_phcorr( :, i, j ), phc0, phc1 ] = ACME( spec_ap_zf_phcorr( :, i, j ), [ 30, 1E-6 ] );
        phc0arr( i, j ) = phc0;
        phc1arr( i, j ) = phc1;
    end
end

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'ACME phase correction applied ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Back to time domain for frequency alignment.                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_ap_zf_phcorr = ifft( ifftshift( spec_ap_zf_phcorr, 1 ) );

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Data returned to time-domain ...');

%% Plot all 2D-CSI data for inspection

count = 1;

hf = figure('Name', 'CSI Matrix Plot', 'MenuBar', 'none', 'NumberTitle', 'off');
ha = axes(hf);
axis off, grid off, box off, hold on

[ Axes, Positions ] = tight_subplot(8, 8, 0.005, 0.005, 0.005);

fmin = floor(min(f_axis_ppm));
fmax = ceil(max(f_axis_ppm));

fmin = 5.0*floor(fmin/5.0);
fmax = 5.0*ceil(fmax/5.0);

ymin = min( real( spec_ap_zf_phcorr( : ) ) );
ymax = max( real( spec_ap_zf_phcorr( : ) ) );

[ ymin, ymax ] = pft_YLimits(ymin, ymax, 0.5);

rr = zeros(1, NPTS, 'int32');
cc = zeros(1, NPTS, 'int32');

for n = 1:NPTS
  tx = coord{n}(1);
  ty = coord{n}(2);
  rr(n) = ty;
  cc(n) = tx;
end

SingleIndex = 1 + 8*(rr - 1) + (cc - 1);    

for row = 1:8
    for col = 1:8 
        set(hf, 'CurrentAxes', Axes(count));
        axis square, box on, grid off, hold on
        if ismember(count, SingleIndex)
          plot( f_axis_ppm, real( spec_ap_zf_phcorr( :, col, row ) ), '-r', 'LineWidth', 1.0 );
        else
          plot( f_axis_ppm, real( spec_ap_zf_phcorr( :, col, row ) ), '-b', 'LineWidth', 1.0 );
        end
        set(Axes(count), 'XDir', 'reverse');
        set(Axes(count), 'XTick', []);
        set(Axes(count), 'XTickLabel', []);
        set(Axes(count), 'YTick', []);
        set(Axes(count), 'YTickLabel', []);   
        xlim( [ fmin, fmax ] );
        ylim( [ ymin, ymax ] );
        if (count == 1)
          text(0.05, 0.90, 'x = 1, y = 1', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
          text(0.05, 0.75, SubFolder, 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');           
        elseif (count == 2)
          text(0.05, 0.90, 'x = 2, y = 1', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
          text(0.05, 0.75, BRUN, 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
        elseif (count == 3)
          text(0.05, 0.90, 'x = 3, y = 1', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
          text(0.05, 0.75, 'PCr at 1.2 ppm', 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
        else
          text(0.05, 0.90, sprintf('%1d, %1d', col, row), 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
        end
        count = count + 1;
    end
end
    
pft_FormatSquarePrinting(hf, ha);

OutputFileNameStub = sprintf('%s - %s - CSI Matrix Plot', SubFolder, BRUN);

pft_ExportGraphsInSeveralFormats(hf, Away, OutputFileNameStub);

pause(0.25);

delete(Axes);

delete(ha);
delete(hf);
    
PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'CSI matrix plot saved ...');

%% IDENTIFY SEPTAL FID's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Take coordinates entered earlier and use them to select FID's from septal voxels.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_septum = zeros( size( fids_ap_zf_phcorr, 1 ), numel( coord ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loop over a cell array of coordinates and create an array of septal FID's.                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : numel( coord )
    xx = coord{ i }( 1 );
    yy = coord{ i }( 2 );
    cc = xx;
    rr = yy;    
    fids_septum( :, i ) = fids_ap_zf_phcorr( :, cc, rr );     % Row, column convention - transposed by PFT
end

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Septal FID''s identified ...');

%% AUTOMATIC FREQUENCY ALIGNMENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Based on the correlation method of Wiegers et al. (MAGMA 2017).                                                                  %
%  Initialise arrays for storing corrected data.                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_septum_fcorr = zeros( size( fids_septum ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate the time axis.                                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_axis = 0 : dw : ( sz( 1 ) - 1 ) * dw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Save calculated terms in a cell array.                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delt = cell( size( fids_septum, 2 ), size( fids_septum, 3 ) ); 
delt_sv = [ 0, 0 ]; % Starting values for delta freq. and delta phase.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Use the highest SNR spectrum in the series as a reference.                                                                       %
%  Test correlation of all other spectra versus the reference.                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : size( fids_septum, 2 )
    
        [ Maximum, Index ] = max( max( abs( fids_septum ) ) ); % Index of the FID with the highest SNR
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine frequency and phase shifts by maximising the cor_fun function %
        % (fminsearch where function has negative sign).                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        delt{ i } = fminsearch( @( delt_sv )cor_fun( delt_sv, ...
                                   fids_septum( :, Index ), ...
                                   fids_septum( :, i ), ...
                                   t_axis ), ...
                                   delt_sv );
                               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply frequency and phase correction terms to FID's.                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fids_septum_fcorr( :, i ) = ...
        fids_septum( :, i ) .* ...
        exp( 1i * 2 * pi * t_axis' * delt{ i }( 1 ) ) .* ...
        exp( 1i * 2 * pi * delt{ i }( 2 ) / 360 );

end

spec_septum_fcorr = fftshift( fft( fids_septum_fcorr), 1 );

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Frequency alignment completed ...');

%% AVERAGE SEPTAL SPECTRA

spec_septum_avg = mean( spec_septum_fcorr, 2 ); 

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Septal spectra averaged ...');

%% Plot selected septal voxels, and their average, for inspection

hf = figure('Name', 'Septal Spectra', 'MenuBar', 'none', 'NumberTitle', 'off');
axis on, box on, grid off, hold on

fmin = floor(min(f_axis_ppm));
fmax = ceil(max(f_axis_ppm));

fmin = 5.0*floor(fmin/5.0);
fmax = 5.0*ceil(fmax/5.0);

ymin = min( real( spec_septum_fcorr( : ) ) );
ymax = max( real( spec_septum_fcorr( : ) ) );

[ ymin, ymax ] = pft_YLimits(ymin, ymax, 0.5);

plot( f_axis_ppm, real( spec_septum_fcorr ), 'LineWidth', 1.0 );
set(gca, 'XDir', 'reverse');
xlim( [ fmin, fmax ] );
ylim( [ ymin, ymax ] );
ht = title(sprintf('%s - %s - Septal Spectra', SubFolder, BRUN));
set(ht, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');
xl = xlabel('^3^1P chemical shift [ppm]'); 
set(xl, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');
yl = ylabel('Signal Intensity [AU]'); 
set(yl, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');

pft_FormatLandscapePrinting(gcf, gca);
OutputFileNameStub = sprintf('%s - %s - Septal Spectra', SubFolder, BRUN);
pft_ExportGraphsInSeveralFormats(gcf, Away, OutputFileNameStub);
pause(0.25);
delete(hf);

hf = figure('Name', 'Averaged Septal Spectrum', 'MenuBar', 'none', 'NumberTitle', 'off');
axis on, box on, grid off, hold on

plot( f_axis_ppm, real( spec_septum_avg ), 'LineWidth', 1.0 );
set(gca, 'XDir', 'reverse');
xlim( [ fmin, fmax ] );
ylim( [ ymin, ymax ] );
ht = title(sprintf('%s - %s - Averaged Septal Spectrum', SubFolder, BRUN));
set(ht, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');
xl = xlabel('^3^1P chemical shift [ppm]'); 
set(xl, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');
yl = ylabel('Signal Intensity [AU]'); 
set(yl, 'FontUnits', 'pixels', 'FontSize', 14, 'FontWeight', 'bold');

pft_FormatLandscapePrinting(gcf, gca);
OutputFileNameStub = sprintf('%s - %s - Averaged Septal Spectrum', SubFolder, BRUN);
pft_ExportGraphsInSeveralFormats(gcf, Away, OutputFileNameStub);
pause(0.25);
delete(hf);

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Septal plots saved ...');

%% Prepare data for export to jMRUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FT spectra prior to saving and split into real and imaginary components.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids_real = real( ifft( ifftshift( spec_septum_avg ) ) );
fids_imag = imag( ifft( ifftshift( spec_septum_avg ) ) );     % Donnie - is this correct ? Yes: 20/12/2018.
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real and imaginary components of spectra.                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spec_real = real( spec_septum_avg );
spec_imag = imag( spec_septum_avg );                          % Donnie - is this correct ? Yes: 20/12/2018.

%% WRITE DATA TO .TXT FILE FOR jMRUI IMPORT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  First, generate .txt. file name using participant ID, series number, and date (in YYYY.MM.DD.HH.mm format) from DICOM file name. %
%  There are two naming conventions for .IMA files, depending on how files are exported.                                            %
%  Test for these to ensure our filename comes out right.                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Then, write a TXT-mode file header.                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( fid, '%s\n\n', 'jMRUI Data Textfile', 'Filename: card31P_recon.txt' );
fprintf( fid, '%s\n', ...
                        [ 'PointsInDataset:', ' ', num2str( size( fids_real, 1 ) ) ], 'DatasetsInFile: 1', ...
                        [ 'SamplingInterval: ', num2str( dw * 1000 ) ], ...
                        'ZeroOrderPhase: 0', ...
                        'BeginTime: 0', ...
                        [ 'TransmitterFrequency:', ' ', num2str( info.csa.ImagingFrequency ), 'E6' ], ...
                        'MagneticField: 3E0', ...
                        'TypeOfNucleus: 1E0', ...
                        [ 'Name of Patient:', ' ', fileName( 1 : mr_ind ) ], ...
                        [ 'Date of Experiment:', ' ', date ], ...
                        'Spectrometer: Siemens_Prisma_3T', ...
                        'AdditionalInfo: na' );
                    
fprintf( fid, '%s', 'SignalNames: Septal_31P_Spectrum' );	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Write the MRS data to text file.                                                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Data written to text-mode jMRUI file ...');

%% Concatenate the 3 plots into a single PDF summary file

Output = sprintf('%s - %s - Cardiac CSI Summary.pdf', SubFolder, BRUN);

if (exist(Output, 'file') == 2)
  delete(Output);
end

Input1 = sprintf('%s - %s - CSI Matrix Plot.pdf', SubFolder, BRUN);
Input2 = sprintf('%s - %s - Septal Spectra.pdf', SubFolder, BRUN);
Input3 = sprintf('%s - %s - Averaged Septal Spectrum.pdf', SubFolder, BRUN);

Inputs = { Input1, Input2, Input3 };

append_pdfs(Output, Inputs{:});

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Plotted graphs concatenated ...');

%% Signal completion

cd(Home);

PHASES = PHASES + 1;

waitbar(double(PHASES)/double(NSTEPS), wb, 'Pre-processing complete !');

pause(1.0);

delete(wb);

end

