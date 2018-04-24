function data_filt = itd_filter( spec_data, sd_filter )
% PURPOSE: Reads in data imported from GE 'P-file', performs Gaussian sliding 
% window filtering in the indirect time-domain (exercise domain) and returns 
% filtered data. Based on work by Rowland et al. (MRM and ISMRM 2016)  

% INPUTS: 'peak_amps' - peak amplitudes in double precision; and 'sd_filter' 
% - an integer representing the standard deviation of the Gaussian sliding 
% window function.
% OUTPUTS: 'filtered_data' - the filtered peak amplitudes, double precision.

% AUTHOR: Donnie Cameron
% DATE: 15/07/2016
% LAST UPDATED: 05/04/2018
%===============================================================

%% CONSTRUCT BLURRING WINDOW
% USER INPUT (for now)
wdth_pts = 9;

% Generate window
window_width = int16( wdth_pts );
half_width = window_width / 2;
alpha = wdth_pts / ( 2 * sd_filter );
gauss_filter = gausswin( wdth_pts, alpha );
gauss_filter = gauss_filter / sum( gauss_filter ); % Normalise.
pad_val = 8;
%wvtool(gausswin(wdthPts,alpha));

%% SMOOTH DATA IN THE INDIRECT TIME-DOMAIN
pad_data = cell( size( spec_data, 1 ), 1 );
pad_data_filt = cell( size( spec_data, 1 ), 1 );
for i = 1 : size( spec_data, 1 )
    pad_data{ i } = padarray( spec_data( i, : ), [ 0 pad_val ], 'symmetric' );
    pad_data_filt{ i } = conv( pad_data{ i }, gauss_filter' );
end

% Cut (transposed) vector down to original size.
data_filt = zeros( size( spec_data ) );
for i = 1 : size( spec_data, 1 ) 
    data_filt( i, : ) = pad_data_filt{ i }( ... 
        pad_val + half_width : end - pad_val - half_width + 1 )';
end

%% PLOT RESULTING SPECTRA
% Plot surface plot of filtered spectra
% figure( 'Name', 'Surface Plot of Filtered 31P Spectra Over Time', ...
%     'NumberTitle', 'off' );
% mesh( real( data_filt ) )
% colormap(autumn)

end
