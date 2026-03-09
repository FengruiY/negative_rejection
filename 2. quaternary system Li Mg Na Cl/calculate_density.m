% Calculate the density of pure water at temp (K)
% Range of validity : [-30 ; 150] ℃

function y = calculate_density(temp)
T = temp-273.15; % ℃

a= -2.8054253*10^-10;
b = 1.0556302*10^-7;
c = -4.6170461*10^-5;
d = -0.0079870401;
e = 16.945176;
f = 999.83952;
g = 0.01687985;

y = (((((a*T+b)*T+c)*T+d)*T+e)*T+f) / (1+g*T); % kg/m3
end