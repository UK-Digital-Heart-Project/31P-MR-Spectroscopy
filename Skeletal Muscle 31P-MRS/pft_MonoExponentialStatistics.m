function S = pft_MonoExponentialStatistics(p)

SteadyStateValue = p(2);
Depletion        = 100.0*p(1)/p(2);
Rate             = p(3);

TimeConstant = 1.0/Rate;

S = { SteadyStateValue, ...
      Depletion, ...
      Rate, TimeConstant }; 
  
end

