function y = solve_kc_alt(kc_conv, jv)
phi = jv*kc_conv.^-1; % Dimensionless

Xi = phi+(1+0.26*phi.^1.4).^-1.7; % Dimensionless

y = kc_conv.*Xi; % m/s
end