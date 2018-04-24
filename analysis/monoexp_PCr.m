function F = monoexpPCr(p,tArr)
% Function for fitting a monoexponential to experimental data in the least-
% squares sense.  
F = p(2)-p(1).*(exp(-p(3).*tArr));
