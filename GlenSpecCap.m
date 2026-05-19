function Result = GlenSpecCap(h)
% Function for calculating a vector of specific water capacities [1/cm]
% from a vector of pressure heads [cm]. The specific capacity function was
% obtained by differentiating the water retention function in Eq. [2] and
% [3] of van Genuchten (1980).
alpha = 0.0115; % [1/cm]
ThetaSat = 0.520;
ThetaRes = 0.218;
n = 2.03;
m = 1-1/(n);
ah = abs(alpha*h);
Result = n*m*alpha*(ThetaSat-ThetaRes)*ah.^(n-1)./(1 + ah.^n).^(m+1);