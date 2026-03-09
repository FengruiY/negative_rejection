function y = solve_kc(N_leaf,W_module,Q0,WR,temp)
Difu_298 = [1.03E-09 7.07E-10 1.33E-09 2.03E-09]';

vis_u = exp(-10.5+530/(temp-146));
Difu = Difu_298.*0.897E-3./298.15.*temp./vis_u;

K = 0.5;
hb = 7e-4; % m
L_mix = 0.006; % m

A = W_module*2*hb*N_leaf; % m2
ub = Q0*(1-WR)/A; % m/s
Pe = 2*hb*ub./Difu; % Dimensionless, vector

rho = calculate_density(temp); % kg/m3
Sc = vis_u/rho./Difu; % Dimensionless, vector

y = 0.753*(K/(2-K))^0.5*Difu/hb.*Sc.^(-1/6).*(Pe*hb/L_mix).^0.5; % mass transfer coefficient k, m/s
end