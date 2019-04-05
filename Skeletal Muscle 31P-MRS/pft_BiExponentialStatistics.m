function S = pft_BiExponentialStatistics(p)

SteadyStateValue = p(2);
Depletion        = 100.0*p(1)/p(2);
Rate1            = p(4);
Rate2            = p(5);

if (Rate1 < Rate2)
  SlowRate       = Rate1;
  FastRate       = Rate2;
  SlowFraction   = 100.0*p(3);
  FastFraction   = 100.0*(1.0 - p(3));
else
  SlowRate       = Rate2;
  FastRate       = Rate1;
  SlowFraction   = 100.0*(1.0 - p(3));
  FastFraction   = 100.0*p(3);
end

SlowTimeConstant = 1.0/SlowRate;
FastTimeConstant = 1.0/FastRate;

S = { SteadyStateValue, ...
      Depletion, ...
      SlowFraction, SlowRate, SlowTimeConstant, ...
      FastFraction, FastRate, FastTimeConstant }; 
  
end

