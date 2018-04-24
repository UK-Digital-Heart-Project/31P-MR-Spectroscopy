function [ results, n_dyn ] = results_parse( file_name )
% PURPOSE: This function finds relevant results from an AMARES results file,
% including PCr and Pi amplitudes, and noise, pH, and [Mg2+] estimates.

% INPUTS: 'file_name' - a string referring to an AMARES results .txt file in 
% the working directory.
% OUTPUTS: 'results' - a table array containing PCr and Pi amplitudes, and 
% noise, pH, and [Mg2+] estimates; and 'n_dyn' - an integer indicating the 
% number of spectra in the results file.

% AUTHOR: Donnie Cameron
% DATE: 05/04/2018
% LAST UPDATED: 05/04/2018

% Check working directory for named file.
d = dir( file_name );
p31_data = fopen( d.name, 'rt' );
% Ensure enough columns for all measured peaks. 
pD = textscan( p31_data, '%s%s%s%s%s%s%s%s%s%s', 'Delimiter', { '\t' } ); 
% WARNING: different no. of fitted metabolites will cause error!
fclose( p31_data );

% Determine number of dynamic spectra by comparing string indices.
I = strcmp( 'Amplitudes (-)', pD{ 1, 1 } );
ind_amp  = find( I );
I = strcmp( 'Standard deviation of Amplitudes (-)', pD{ 1, 1 } );
ind_amp_sd  = find( I );
% No. of spectra = number of lines between 'Amp' and 'SD Amp' headings - 1. 
n_dyn = ind_amp_sd - ind_amp - 1;

% Store PCr and Pi peak amplitudes and SD of PCr amplitude. 
pcr_amp = str2double( pD{ 1, 2 }( ind_amp + 1 : ind_amp + n_dyn ) );
pi_amp = str2double( pD{ 1, 1 }( ind_amp + 1 : ind_amp + n_dyn ) );
pcr_amp_sd = str2double( pD{ 1, 2 }( ind_amp_sd + 1 : ind_amp_sd + n_dyn ) );

% Store noise estimates for SNR calculation later.
I = strcmp( 'Noise : ', pD{ 1, 1 } );
ind_noise  = find( I );
noise_dyn = str2double( pD{ 1, 1 }( ind_noise + 1 : ind_noise + n_dyn ) );

% Store estimated pH values to check for acidosis later.
I = strcmp( 'pH', pD{ 1, 1 } );
ind_pH  = find( I );
pH_dyn = str2double( pD{ 1, 1 }( ind_pH + 1 : ind_pH + n_dyn ) );

% Store Mg2+ concentrations for use in future work.
I = strcmp( '[Mg2+]', pD{ 1, 3 } );
ind_mg  = find( I );
mg_dyn = str2double( pD{ 1, 1 }( ind_mg + 1 : ind_mg + n_dyn ) );

% Save results as a table array.
results = table( pcr_amp, ...
    pi_amp, ...
    pcr_amp_sd, ...
    noise_dyn, ...
    pH_dyn, ...
    mg_dyn );

end
