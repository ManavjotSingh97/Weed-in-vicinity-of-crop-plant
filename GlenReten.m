function Result = GlenReten(h)
% Function for calculating a vector of soil volumetric water contents
% [cm^3/cm^3] from a vector of pressure heads [cm]. The water retention
% function is from Eq. [2] and [3] of van Genuchten (1980).
alpha = 0.0115; % [1/cm]
ThetaSat = 0.520;
ThetaRes = 0.218;
n = 2.03;
m = 1-1/(n);
Result = ThetaRes + ((ThetaSat-ThetaRes)./((1 + abs(alpha*h).^n).^m));