function F = biexpPCr(p,tArr)
% Function for fitting a biexponential to experimental data in the least-
% squares sense.  
F = p(2)-p(1).*(p(3)*exp(-p(4).*tArr)+(1-p(3)).*exp(-p(5).*tArr));
