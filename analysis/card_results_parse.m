function results = card_results_parse( file_name )
%% PURPOSE: This function finds relevant results from an AMARES results file,
% including PCr, y-ATP, and DPG amplitudes, and noise.

% INPUTS: 'file_name' - a string referring to an AMARES results .txt file in 
% the working directory.
% OUTPUTS: 'results' - a table array containing PCr, y-ATP, & DPG amplitudes,
% and noise estimates.

% AUTHOR: Donnie Cameron
% DATE: 05/04/2018
% LAST UPDATED: 03/06/2018

% Check working directory for named file.
d = dir( file_name );
p31_data = fopen( d.name, 'rt' );
% Ensure enough columns for all measured peaks. 
pD = textscan( p31_data, '%s%s%s%s%s%s%s%s%s%s%s', 'Delimiter', { '\t' } ); 
% WARNING: different no. of fitted metabolites will cause error!
fclose( p31_data );

% Determine number of dynamic spectra by comparing string indices.
I = strcmp( 'Amplitudes (-)', pD{ 1, 1 } );
ind_amp  = find( I );
I = strcmp( 'Standard deviation of Amplitudes (-)', pD{ 1, 1 } );
ind_amp_sd  = find( I );

% Store PCr, y-ATP, and DPG peak amplitudes and their SDs. 
pcr_amp = str2double( pD{ 1, 1 }( ind_amp + 1 ) );
pcr_amp_sd = str2double( pD{ 1, 1 }( ind_amp_sd + 1 ) );
yatp1_amp = str2double( pD{ 1, 2 }( ind_amp + 1 ) );
yatp1_amp_sd = str2double( pD{ 1, 2 }( ind_amp_sd + 1 ) );
yatp2_amp = str2double( pD{ 1, 3 }( ind_amp + 1 ) );
yatp2_amp_sd = str2double( pD{ 1, 3 }( ind_amp_sd + 1 ) );
dpg2_amp = str2double( pD{ 1, 9 }( ind_amp + 1 ) );
dpg2_amp_sd = str2double( pD{ 1, 9 }( ind_amp_sd + 1 ) );
dpg3_amp = str2double( pD{ 1, 10 }( ind_amp + 1 ) );
dpg3_amp_sd = str2double( pD{ 1, 10 }( ind_amp_sd + 1 ) );

% Store noise estimates for SNR calculation later.
I = strcmp( 'Noise : ', pD{ 1, 1 } );
ind_noise  = find( I );
noise = str2double( pD{ 1, 1 }( ind_noise + 1 ) );

% Save results as a table array.
results = table( pcr_amp, ...
    yatp1_amp, ...
    yatp2_amp, ...
    dpg2_amp, ...
    dpg3_amp, ...
    pcr_amp_sd, ...
    yatp1_amp_sd, ...
    yatp2_amp_sd, ...
    dpg2_amp_sd, ...
    dpg3_amp_sd, ...
    noise );

end