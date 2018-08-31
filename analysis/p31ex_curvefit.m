function [ results_table, fits ] = p31ex_curvefit( pcr_amp, pcr_init, ...
    pcr_min_index, t_arr_init, t_dyn, bout_ind )
% PURPOSE: This function fits 31P Phosphocreatine time course data for 
% estimation of kPCr.

% INPUTS: 'pcr_amp' - vector containing PCr peak amplitudes; 
%   'pcr_init' - vector containing PCr amplitudes pre-normalisation;
%   'pcr_min_index' - index of PCr minimum, assumed to be exercise end pt.;
%   't_arr_init' - vector, the time axis of the data; and 
%   't_dyn' - repetition time of MRS sequence * # of excitations per dynamic.
% OUTPUTS: 'results_table' - table array containing the estimated kPCr values 
% from monoexponential and biexponential fits; and
%   'fits' - vector of strings indicating the best fit for kPCr from mono-
% exponential and biexponential models and the index of the fitting start pt.

% AUTHOR: Donnie Cameron
% DATE: 06/04/2018
% LAST UPDATED: 31/08/2018
%==============================================

% SET STOP PT FOR ITERATIVE CURVE-FITTING
stop_pt = pcr_min_index + 2;  % Search up to ex end pt plus 2

% Don't display curve-fitting messages - speeds up fitting.
options = optimset( 'display', 'off' ); 

% Initial starting values and bounds for mono- and biexponential fits.
sv_mono = [ 0, 0, 0.03 ];
lb_mono = [ 0, 0, 0 ];
ub_mono = [ inf, inf, inf ];

sv_bi = [ 0, 0, 0.5, 0.03, 0 ];
lb_bi = [ 0, 0, 0, 0, 0 ];
ub_bi = [ inf, inf, 1, inf, inf ];

% Loop over all elements in PCr array, selecting different start pts. for
% the curve fit
alpha = max( pcr_init( 1 : 10 ) ); % Alpha represents 'PCr rest', which can be 
                                   % estimated from pre-/post- exercise data.
if max( pcr_init( end - 9 : end ) ) > alpha
    alpha = max( pcr_init( end - 9 : end ) ); % If PCr val. > PCr rest after
end                                           % exercise, calculate new
                                              % alpha from last 10 pts. 

% Preallocate prior to 'for' loop - for speed.                                               
t_arr = cell( stop_pt, 1 );

param_mono = cell( stop_pt, 1 );
res_norm_mono = zeros( stop_pt, 1 );
RMSE_mono = zeros( stop_pt, 1 );
kPCr_mono = zeros( stop_pt, 1 );

param_bi = cell( stop_pt, 1 );
resnorm_bi = zeros( stop_pt, 1 );
RMSE_bi = zeros( stop_pt, 1 );
kPCr_bi1 = zeros( stop_pt, 1 );
kPCr_bi2 = zeros( stop_pt, 1 );

% Loop over start values on the recovery curve. Don't fit points on plateau.                                                
for i = 1 : stop_pt 
    pcr = pcr_amp( i : end );
    t_arr{ i } = ones( numel( pcr ), 1 );  % Store timings for init. experiment.
    t_arr{ i } = find( t_arr{ i } ) .* t_dyn;
    t_arr2 = t_arr{ i }; % Avoids errors when referencing cells in lsqcurvefit.    
    beta = min( pcr );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Fit 'monoexponential rise to a maximum' to the data. Test for pts
    % adjacent to the minimum.
    sv_mono( 1 ) = beta;
    sv_mono( 2 ) = alpha;
    % Least-squares curve_fitting
    [ param_mono{ i }, res_norm_mono( i ) ] = lsqcurvefit( ...
        @( svMono, tArr2 )monoexp_pcr( svMono, tArr2 ), ...
        sv_mono, ...
        t_arr2, pcr,...
        lb_mono, ub_mono, ...
        options );
    % Calculate root mean square error of fit.
    RMSE_mono( i ) = sqrt( res_norm_mono( i ) / ... 
        ( numel( t_arr2 ) - numel( param_mono{ i } ) ) );
    
    % Store monoexp. kPCr estimate for easy plotting later.
    kPCr_mono( i ) = param_mono{ i }( 3 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit biexponential to the data. Again, test for pts adjacent to minimum.
    sv_bi( 1 ) = beta;
    sv_bi( 2 ) = alpha;
    % Least-squares curve_fitting
    [ param_bi{ i }, resnorm_bi( i ) ] = lsqcurvefit( ...
        @( svBi, tArr2 )biexp_pcr( svBi, tArr2 ), ...
        sv_bi, ...
        t_arr2, pcr, ...
        lb_bi, ub_bi, ...
        options );
    % Calculate root mean square error of fit.
    RMSE_bi( i ) = sqrt( resnorm_bi( i ) / ...
        ( numel( t_arr2 ) - numel( param_bi{ i } ) ) );
    
    % Store biexp. kPCr1 estimate for easy plotting later.
    kPCr_bi1( i ) = param_bi{ i }( 4 ); 
    % Store biexp. kPCr2 estimate for easy plotting later.
    kPCr_bi2( i ) = param_bi{ i }( 5 );
    
end

% Pick best fit based on RMSE
[ ~, J ] = min( RMSE_mono );
% Pick best fit based on RMSE
[ ~, K ] = min( RMSE_bi ); 

%% Plot best fits for both models.
h = figure( 'Name', [ 'Mono- and Biexponential Fit Comparison, Bout ', ...
    num2str( bout_ind ) ], 'NumberTitle', 'off' );
F1 = biexp_pcr( param_bi{ K }, t_arr{ K } );
F2 = monoexp_pcr( param_mono{ J }, t_arr{ J } );
scatter( t_arr_init, pcr_init, 'k', '.' );
hold on
fitPlots = plot( t_arr{ K } + ( K - 1 ) .* t_dyn, F1, ...
    t_arr{ J } + ( J - 1 ) .* t_dyn, F2 );
fitPlots( 1 ).LineWidth = 2;
fitPlots( 2 ).LineWidth = 2;
hold off
xlabel( 'Time (s)' )
ylabel( 'PCr Signal Intensity (a.u.)' )
legend( 'Signal', 'Biexponential Fit', 'Monoexponential Fit', ...
    'Location', 'southeast' )
xlim( [0 max( t_arr_init ) ] );
ylim( [ 0.95 * min( pcr_init ) 1.02 * max( pcr_init ) ] ); 

% Time constant of monoexponential fit.
mono_KPCr = param_mono{ J }( 3 );
% 'Fast' and 'slow' biexponential parameters sometimes swap. Pick the one
% that most resembles the monoexponential time constant.
bi_tc = [ param_bi{ K }( 4 ) param_bi{ K }( 5 ) ];
[ ~, I ] = min( abs( bi_tc - mono_KPCr ) );
bi_KPCr = bi_tc( I );
fits{ 1 } = sprintf( 'Monoexponential kPCr = %s at dynamic %s (%s s)\n', ...
    num2str( mono_KPCr ), num2str( J ), num2str( J * 3 ) );
fits{ 2 } = sprintf( 'Biexponential kPCr(1) = %s at dynamic %s (%s s)\n', ...
    num2str( bi_KPCr ), num2str( K ), num2str( J * 3 ) );

%% Create results table.
results_table = table( mono_KPCr, ...
    bi_KPCr );
     
end
