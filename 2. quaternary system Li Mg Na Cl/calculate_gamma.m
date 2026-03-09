function gamma = calculate_gamma(c,epsilon,temp)
z = [1 2 1 -1]';

c = c./1000;
I = 0;
for i = 1:length(c)
    I = I+0.5*z(i)^2*c(i);
end

NA = 6.022e23;
e = 1.602e-19;
epsilon_0 = 8.854e-12;
kB = 1.380649e-23;
rho = calculate_density(temp);

A = rho^0.5*((2*pi*NA)^0.5/log(10))*(e^2/(4*pi*epsilon_0*epsilon*kB*temp))^1.5;

gamma = 10.^(-A*(sqrt(I)/(1+sqrt(I))-0.3*I)*z.^2);
end