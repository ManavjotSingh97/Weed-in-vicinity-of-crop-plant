function Result = GlenCond(h)
% Function for calculating a vector of hydraulic conductivities [cm/hr]
% from a vector of pressure heads [cm]. The conductivity function is from
% Eq. [9] of van Genuchten (1980).
alpha = 0.0115; % [1/cm]
Ksat = 31.6/24; % [cm/hr]
n = 2.03;
m = 1-1/(n);
ah = abs(alpha*h);
if h < 0
    Result = Ksat * (1 - ah.^(n-1).*(1+ah.^n).^(-m)).^2./(1+ah.^n).^(m/2);
else 
    Result = Ksat;
end