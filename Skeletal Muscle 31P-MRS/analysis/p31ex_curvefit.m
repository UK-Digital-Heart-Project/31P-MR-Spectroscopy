function [ MonoExponentialStatistics, BiExponentialStatistics, DCStatistics, MonoKPCr, BiKPCr ] = ...
p31ex_curvefit( pcr_amp, pcr_init, pcr_min_index, t_arr_init, t_dyn, bout_ind, Title )
% PURPOSE: This function fits 31P Phosphocreatine time course data for 
% estimation of kPCr.

% INPUTS: 'pcr_amp' - vector containing PCr peak amplitudes; 
%   'pcr_init' - vector containing PCr amplitudes pre-normalisation;
%   'pcr_min_index' - index of PCr minimum, assumed to be exercise end pt.;
%   't_arr_init' - vector, the time axis of the data; and 
%   't_dyn' - repetition time of MRS sequence * # of excitations per dynamic.

% AUTHOR: Donnie Cameron
% DATE: 06/04/2018
% LAST UPDATED: 31/08/2018
%==============================================

% SET STOP PT FOR ITERATIVE CURVE-FITTING
stop_pt = pcr_min_index + 2;  % Search up to ex end pt plus 2

% Don't display curve-fitting messages - speeds up fitting.
options = optimset( 'display', 'off' ); 

% Loop over all elements in PCr array, selecting different start pts. for the curve fit
alpha = max( pcr_init( 1 : 10 ) ); % Alpha represents 'PCr rest', which can be estimated from pre-/post- exercise data.
if max( pcr_init( end - 9 : end ) ) > alpha
    alpha = max( pcr_init( end - 9 : end ) ); % If PCr val. > PCr rest after exercise, calculate new alpha from last 10 pts.
end      

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
    % Fit 'monoexponential rise to a maximum' to the data. Test for pts adjacent to the minimum.
    % sv_mono( 1 ) = beta;
    % sv_mono( 2 ) = alpha;
    
    sv_mono = [ max(pcr) - min(pcr), max(pcr), 0.03 ];
    lb_mono = [ 0.0, 0.0, 0.0 ];
    ub_mono = [ 0.5*max(pcr), 2.0*max(pcr), 1.0 ];
    
    % Least-squares curve_fitting
    [ param_mono{ i }, res_norm_mono( i ) ] = lsqcurvefit( @( svMono, tArr2 )monoexp_PCr( svMono, tArr2 ), sv_mono, t_arr2, pcr, lb_mono, ub_mono, options );
    % Calculate root mean square error of fit.
    RMSE_mono( i ) = sqrt( res_norm_mono( i ) / ( numel( t_arr2 ) - numel( param_mono{ i } ) ) );
    
    % Store monoexp. kPCr estimate for easy plotting later.
    kPCr_mono( i ) = param_mono{ i }( 3 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit biexponential to the data. Again, test for pts adjacent to minimum.
    % sv_bi( 1 ) = beta;
    % sv_bi( 2 ) = alpha;
    sv_bi = [ max(pcr) - min(pcr), max(pcr), 0.5, 0.02, 0.05 ];
    lb_bi = [ 0.0, 0.0, 0.0, 0.0, 0.0 ];
    ub_bi = [ 0.5*max(pcr), 2.0*max(pcr), 1.0, 1.0, 1.0 ];    
    
    % Least-squares curve_fitting
    [ param_bi{ i }, resnorm_bi( i ) ] = lsqcurvefit( @( svBi, tArr2 )biexp_PCr( svBi, tArr2 ), sv_bi, t_arr2, pcr, lb_bi, ub_bi, options );
    % Calculate root mean square error of fit.
    RMSE_bi( i ) = sqrt( resnorm_bi( i ) / ( numel( t_arr2 ) - numel( param_bi{ i } ) ) );
    
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
h = figure( 'Name', [ 'Mono-Exponential and Bi-Exponential Fit Comparison, Bout ', num2str( bout_ind ) ], 'NumberTitle', 'off', 'MenuBar', 'none' );
axis on, box on, hold on
F1 = biexp_PCr( param_bi{ K }, t_arr{ K } );
F2 = monoexp_PCr( param_mono{ J }, t_arr{ J } );
scatter( t_arr_init, pcr_init, 25, 'k', 'filled' );
hold on
fitPlots = plot( t_arr{ K } + ( K - 1 ) .* t_dyn, F1, t_arr{ J } + ( J - 1 ) .* t_dyn, F2 );
fitPlots( 1 ).LineWidth = 2;
fitPlots( 2 ).LineWidth = 2;
fitPlots( 1 ).Color = 'g';
fitPlots( 2 ).Color = 'm';
hold off
xlabel( 'Time / sec' )
ylabel( 'PCr Signal Intensity / A.U.' )
title(Title);
legend( 'Signal', 'Bi-Exponential Fit', 'Mono-Exponential Fit', 'Location', 'SouthEast' )
legend boxoff
xlim([0, 50.0*ceil(max(t_arr_init)/50.0)]);
[ MinY, MaxY ] = pft_YLimits(min(pcr_init), max(pcr_init), 0.2);
ylim([MinY, MaxY]);

pft_FormatLandscapePrinting(gcf, gca);

pft_ExportGraphsInSeveralFormats(gcf, pwd, sprintf('Bout %1d Curve-Fitting', bout_ind));

delete(h);

% Time constant of monoexponential fit.
mono_KPCr = param_mono{ J }( 3 );

% 'Fast' and 'slow' biexponential parameters sometimes swap, so pick the one that most resembles the monoexponential time constant
bi_tc = [ param_bi{ K }( 4 ) param_bi{ K }( 5 ) ];
[ ~, I ] = min( abs( bi_tc - mono_KPCr ) );
bi_KPCr = bi_tc( I );

%% Create the summary statistics for the two fits
MonoExponentialStatistics = pft_MonoExponentialStatistics(param_mono{ J });

BiExponentialStatistics = pft_BiExponentialStatistics(param_bi{ K });

DCStatistics = { mono_KPCr, J, 3*J, bi_KPCr, K, 3*K };

%% Create results table and pass out some standalone copies of the rate constants

MonoKPCr = mono_KPCr;
BiKPCr   = bi_KPCr;
     
end
