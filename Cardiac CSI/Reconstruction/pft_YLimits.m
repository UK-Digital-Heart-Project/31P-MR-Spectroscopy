function [ Mini, Maxi ] = pft_YLimits(LoY, HiY, Factor)

Maxi = HiY;
Step = Factor*10^floor(log10(Maxi));
Maxi = Step*ceil(Maxi/Step);

Mini = LoY;
Mini = Step*floor(Mini/Step);

end

