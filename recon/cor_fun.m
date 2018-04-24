function F = cor_fun( delt, fid_ref, fid_n, t_axis )
% PURPOSE: Frequency and phase alignment of an FID with a reference FID by 
% maximising % the correlation between the two in the frequency domain. 
% Based on work by Wiegers et al. (MAGMA 2017)  

% INPUTS: 'delt' (vector) holds the minimisation terms, where delt( 1 ) is the 
% frequency shift in Hz and delt( 2 ) is the phase shift in degrees; 
%   'fid_ref' (complex vector) is the reference FID that the other FID will 
% be aligned to; 
%   'fid_n' (complex vector) is the FID that is to be aligned to 'fid_ref'; &
%   't-axis' (vector) is the time axis of the time-domain data.
% OUTPUTS: 'F' (scalar) is the result of the minimisation.

% AUTHOR: Donnie Cameron
% DATE: 09/04/2018
% LAST UPDATED: 10/04/2018
%===============================================================
fid_n = fid_n .* exp( 1i * 2 * pi * t_axis' * delt( 1 ) ) * ...
    exp( 1i * 2 * pi * delt( 2 ) / 360 );

spec_ref = fft( fid_ref );
spec_n = fft( fid_n );

F = - real( dot( spec_n, spec_ref ) ) / ( norm( spec_ref ) * norm( spec_n ) );
