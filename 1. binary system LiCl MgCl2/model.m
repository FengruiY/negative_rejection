% Function to carry out coupon-scale DSPM-DE simulation based on: 
% Geraldes, V. & Brites Alves, A. M. Computer program for simulation of mass transport in nanofiltration membranes. 
% Journal of Membrane Science 321, 172–182 (2008).

% membrane_paras = [rp, xe, cx, epsilon_p]'
% operation_paras = {jv, temp, cb}'
% other_paras = {kc, eta, tolerance, N}'

function y = model(membrane_paras,operation_paras,other_paras)
%% Read parameters
rp = membrane_paras(1);
xe = membrane_paras(2);
cx = membrane_paras(3);
epsilon_p = membrane_paras(4);

jv = operation_paras{1};
temp = operation_paras{2};
cb = operation_paras{3};
n = length(cb);

kc = other_paras{1};
eta = other_paras{2};
tolerance = other_paras{3};
N = other_paras{4};

%% Define other parameters and constants
% Constants
F = 96485.3321233100184; % C/mol
R = 8.3145; % J/(K·mol)
e0 = 1.602e-19; % electronic charge, C
epsilon_0 = 8.854e-12; % vacuum permittivity, F/m
kB = 1.380649e-23; % Boltzmann constant, J/K

% Other parameters
z = [1 2 -1]';
Difu_298 = [1.03E-09 7.07E-10 2.03E-09]'; % Ion diffusivity at 25℃, m2/s

% Calibrating ion diffusivity at given temperature
vis_u = exp(-10.5+530/(temp-146)); % μ dynamic viscosity, Pa·s
Difu = Difu_298.*0.897E-3./298.15.*temp./vis_u;

r = [2.38E-10 3.50E-10 1.21E-10]'; % Stokes radii of ions, m
lambda = r./rp;

Phi = zeros(n,1); % Steric partitioning factor Φ_i
for i = 1:n
    if lambda(i) >= 1
        Phi(i) = 0;
    else
        Phi(i) = (1-lambda(i))^2;
    end
end

epsilon_b = 78.4;
delta_W = (z.^2./r.*e0^2./(8*pi*epsilon_0)).*(1/epsilon_p-1/epsilon_b);

PhiB = exp(-delta_W./(kB*temp)); % Dielectric exclusion factor, Φ_B,i

% ----------------------------------- Starting calculation below -----------------------------------
%% Initialising equation matrices
% x = zeros((n+1)*(N+2),1);
b = zeros((n+1)*(N+2),1);

x_old = zeros((n+1)*(N+2),1);

rows = zeros(n*(6*N+5),1);
cols = zeros(n*(6*N+5),1);
vals = zeros(n*(6*N+5),1);

%% Initialising the first guess of variables
for i=1:n
    x_old(1+(i-1)*(N+2):i*(N+2)) = cb(i);
end

x_old(n*(N+2)+1:end,1) = 0;

%% Main loop
convergence = false;
ii = 1; % Step number indicator
while convergence == false
    if ii > 50000
        error('Failed to converge, reduce eta.');
    end

    %% Using eq. 22 to update A and b
    count = 1;
    for i = 1:n
        rows(count) = i;
        cols(count) = 1+(i-1)*(N+2);
        vals(count) = kc(i)-jv;
        count = count+1;

        rows(count) = i;
        cols(count) = i*(N+2);
        vals(count) = jv;
        count = count+1;

        rows(count) = i;
        cols(count) = (n+1)*(N+2);
        vals(count) = z(i)*x_old(1+(i-1)*(N+2))*Difu(i)*F/(R*temp);
        count = count+1;

        b(i) = kc(i)*cb(i);
    end

    %% Using eq. 23 to update A and b
    for i = 1:n
        rows(count) = n+1;
        cols(count) = 1+(i-1)*(N+2);
        vals(count) = z(i);
        count = count+1;
    end

    % b(n+1) = 0;

    %% Using eq. 24 to update A and b
    c_i1 = x_old(2+(0:n-1)*(N+2));
    gamma_i1 = calculate_gamma(c_i1,epsilon_p,temp);

    c_im = x_old(1+(0:n-1)*(N+2));
    gamma_im = calculate_gamma(c_im,epsilon_b,temp);

    for i = 1:n
        rows(count) = n+1+i;
        cols(count) = 2+(i-1)*(N+2);
        vals(count) = gamma_i1(i);
        count = count+1;

        rows(count) = n+1+i;
        cols(count) = n*(N+2)+1;
        vals(count) = gamma_im(i)*x_old(1+(i-1)*(N+2))*(z(i)*F/(R*temp))*Phi(i)*PhiB(i)*exp(-z(i)*F*x_old(n*(N+2)+1)/(R*temp));
        count = count+1;

        b(n+1+i) = gamma_im(i)*x_old(1+(i-1)*(N+2))*Phi(i)*PhiB(i)*exp(-z(i)*F*x_old(n*(N+2)+1)/(R*temp))*(1+z(i)*F*x_old(n*(N+2)+1)/(R*temp));
    end

    %% Using eq. 25 to update A and b
    [Di_p,Ki_c] = solve_K_D(lambda,Difu,n);
    delta_xj = xe/(N-1);

    for i = 1:n
        for j = 1:N-1
            % A
            rows(count) = 2*n+1+(i-1)*(N-1)+j;
            cols(count) = 1+j+(i-1)*(N+2);
            vals(count) = Di_p(i)/delta_xj+0.5*Ki_c(i)*jv-0.5*z(i)*Di_p(i)*(F/(R*temp))* ((x_old(n*(N+2)+j+1)-x_old(n*(N+2)+j))/delta_xj) ;
            count = count+1;

            % B
            rows(count) = 2*n+1+(i-1)*(N-1)+j;
            cols(count) = 2+j+(i-1)*(N+2);
            vals(count) = -Di_p(i)/delta_xj+0.5*Ki_c(i)*jv-0.5*z(i)*Di_p(i)*(F/(R*temp))* ((x_old(n*(N+2)+j+1)-x_old(n*(N+2)+j))/delta_xj) ;
            count = count+1;

            % C
            rows(count) = 2*n+1+(i-1)*(N-1)+j;
            cols(count) = n*(N+2)+j;
            vals(count) = 0.5*z(i)*Di_p(i)*(F/(R*temp))* ((x_old((i-1)*(N+2)+1+j+1)+x_old((i-1)*(N+2)+1+j))/delta_xj) ;
            count = count+1;

            % D
            rows(count) = 2*n+1+(i-1)*(N-1)+j;
            cols(count) = n*(N+2)+j+1;
            vals(count) = -0.5*z(i)*Di_p(i)*(F/(R*temp))* ((x_old((i-1)*(N+2)+1+j+1)+x_old((i-1)*(N+2)+1+j))/delta_xj) ;
            count = count+1;

            % E
            rows(count) = 2*n+1+(i-1)*(N-1)+j;
            cols(count) = i*(N+2);
            vals(count) = -jv;
            count = count+1;

            % F
            b(2*n+1+(i-1)*(N-1)+j) = -0.5*z(i)*Di_p(i)*(F/(R*temp))* (x_old((i-1)*(N+2)+1+j+1)+x_old((i-1)*(N+2)+1+j))* ((x_old(n*(N+2)+j+1)-x_old(n*(N+2)+j))/delta_xj) ;
        end
    end

    %% Using eq. 32 to update A and b
    for j=1:N
        for i=1:n
            rows(count) = n+1+n*N+j;
            cols(count) = 1+j+(i-1)*(N+2);
            vals(count) = z(i);
            count = count+1;
        end
        b(n+1+n*N+j) = -cx;
    end

    %% Using eq. 33 to update A and b
    c_iN = x_old(1+N+(0:n-1)*(N+2));
    gamma_iN = calculate_gamma(c_iN,epsilon_p,temp);

    c_ip = x_old((1:n)*(N+2));
    gamma_ip = calculate_gamma(c_ip,epsilon_b,temp);

    for i = 1:n
        rows(count) = n+1+n*N+N+i;
        cols(count) = 1+N+(i-1)*(N+2);
        vals(count) = gamma_iN(i);
        count = count+1;

        rows(count) = n+1+n*N+N+i;
        cols(count) = n*(N+2)+N;
        vals(count) = gamma_ip(i)*x_old(i*(N+2))*z(i)*F/(R*temp)*Phi(i)*PhiB(i)*exp(z(i)*F/(R*temp)*(x_old((n+1)*(N+2)-1)-x_old(n*(N+2)+N)));
        count = count+1;

        rows(count) = n+1+n*N+N+i;
        cols(count) = (n+1)*(N+2)-1;
        vals(count) = -gamma_ip(i)*x_old(i*(N+2))*z(i)*F/(R*temp)*Phi(i)*PhiB(i)*exp(z(i)*F/(R*temp)*(x_old((n+1)*(N+2)-1)-x_old(n*(N+2)+N)));
        count = count+1;

        b(n+1+n*N+N+i) = gamma_ip(i)*x_old(i*(N+2))*Phi(i)*PhiB(i)*exp(z(i)*F/(R*temp)*(x_old((n+1)*(N+2)-1)-x_old(n*(N+2)+N))) * (1-z(i)*F/(R*temp)*(x_old((n+1)*(N+2)-1)-x_old(n*(N+2)+N)));
    end

    %% Using eq. 34 to update A and b
    for i = 1:n
        rows(count) = (n+1)*(N+2);
        cols(count) = i*(N+2);
        vals(count) = z(i);
        count = count+1;
    end

    % b((n+1)*(N+2)) = 0;

    %% Constructing sparse matrix A
    A = sparse(rows,cols,vals,(n+1)*(N+2),(n+1)*(N+2)); %%%
    
    % Define A as a dense matrix:
    % A = zeros((n+1)*(N+2),(n+1)*(N+2));
    % linearIdx = sub2ind([(n+1)*(N+2),(n+1)*(N+2)], rows, cols);
    % A(linearIdx) = vals;

    %% Solve new x
    % Set the warning for determining whether A is close to a singular matrix to an error
    warning('error', 'MATLAB:singularMatrix');
    warning('error', 'MATLAB:nearlySingularMatrix');

    x = A\b;

    %% Check the validity of x
    % Check if the solution contains complex numbers
    if ~isreal(x)
        error('The solution x contains complex numbers.'); % Usually caused by a negative input value for jv
    end

    ii = ii+1;

    %% Convergence check
    % 1. Normalized residual of the discretized Nernst–Planck equation for each ion
    R_NP = zeros(n,1);
    RJT = 0;

    for i = 1:n
        RJT = RJT+abs(jv*x(i*(N+2)));

        for j = 1:N-1
            delta_ci = x((i-1)*(N+2)+1+j+1)-x((i-1)*(N+2)+1+j); % c(i,j+1)-c(i,j)
            add_ci = x((i-1)*(N+2)+1+j+1)+x((i-1)*(N+2)+1+j); % c(i,j+1)+c(i,j)
            delta_phi = x(n*(N+2)+j+1)-x(n*(N+2)+j); % φ(j+1)-φ(j)

            R_NP(i) = R_NP(i)+abs( ...
                -Di_p(i)*delta_ci/delta_xj-0.5*z(i)*add_ci*Di_p(i)*F/(R*temp)*delta_phi/delta_xj+0.5*Ki_c(i)*add_ci*jv-jv*x(i*(N+2)) ...
                );
        end

        if Phi(i) < 1e-15
            R_NP(i) = 0;
        end
    end

    R_NP = R_NP./RJT;

    % 2. Equilibrium equations residuals
    c_i1 = x_old(2+(0:n-1)*(N+2));
    c_i1_new = x(2+(0:n-1)*(N+2));
    gamma_i1 = calculate_gamma(c_i1,epsilon_p,temp);

    c_im = x_old(1+(0:n-1)*(N+2));
    c_im_new = x(1+(0:n-1)*(N+2));
    gamma_im = calculate_gamma(c_im,epsilon_b,temp);

    c_iN = x_old(1+N+(0:n-1)*(N+2));
    c_iN_new = x(1+N+(0:n-1)*(N+2));
    gamma_iN = calculate_gamma(c_iN,epsilon_p,temp);

    c_ip = x_old((1:n)*(N+2));
    c_ip_new = x((1:n)*(N+2));
    gamma_ip = calculate_gamma(c_ip,epsilon_b,temp);

    R_EQ = 0;
    for i = 1:n
        if Phi(i) > 1e-10 && c_im_new(i) > 1e-11
            R_EQ = R_EQ+abs( ...
                c_i1_new(i)*gamma_i1(i)/(c_im_new(i)*gamma_im(i)*Phi(i)*PhiB(i)*exp(-z(i)*F/(R*temp)*x(n*(N+2)+1)))-1 ...
                );
            if c_ip_new(i) > 1e-11
                R_EQ = R_EQ+abs(...
                    c_iN_new(i)*gamma_iN(i)/(c_ip_new(i)*gamma_ip(i)*Phi(i)*PhiB(i)*exp(-z(i)*F/(R*temp)*(x(n*(N+2)+N)-x((n+1)*(N+2)-1))))-1 ...
                    );
            end
        end
    end

    % 3. Feed/membrane mass-transfer equation residuals
    R_MS = 0;

    for i = 1:n
        R_MS = R_MS+abs(-jv*c_ip_new(i)-kc(i)*(c_im_new(i)-cb(i))+jv*c_im_new(i)-z(i)*F/(R*temp)*Difu(i)*c_im_new(i)*x((n+1)*(N+2)));
    end

    R_MS = R_MS/RJT;

    % 4. Electroneutrality residuals
    R_EN1 = 0;
    for j = 1:N
        RR = 0;
        for i = 1:n
            RR = RR+z(i)*x((i-1)*(N+2)+1+j);
        end
        R_EN1 = R_EN1+abs(RR+cx);
    end

    R_EN2 = 0;
    for i=1:n
        R_EN2 = R_EN2+z(i)*c_ip_new(i);
    end
    R_EN2 = abs(R_EN2);

    R_EN_1 = 0;
    for i = 1:n
        R_EN_1 = R_EN_1+abs(z(i)*cb(i));
    end
    R_EN_1 = R_EN_1+abs(cx);

    R_EN = (R_EN1+R_EN2)/R_EN_1;

    % Check convergence
    if all(R_NP < tolerance) && R_EQ < tolerance && R_MS < tolerance && R_EN < tolerance
        convergence = true;
    end

    %% Under-relax the computed solutions and update
    % 1. Update phi
    rangePhi = (n*(N+2)+1) : ((n+1)*(N+2)-1);
    x_old(rangePhi) = x_old(rangePhi)+eta*(x(rangePhi)-x_old(rangePhi));

    % 2. Update ci
    for i = 1:n*(N+2)
        % Check if x_old(i) == 0 or (x(i) - x_old(i)) == 0
        EPS = 1e-14;
        if abs(x(i) - x_old(i)) < EPS || abs(x_old(i)) < EPS
            % Skip the update
            % x_old(i) = x_old(i)
            continue;
        end

        delta_cr = min(abs(x_old(i)/(x(i)-x_old(i))),1)*(x(i)-x_old(i))/x_old(i);
        x_old(i) = x_old(i)*(1+eta*delta_cr);
    end

    % 3. Update xi without under-relaxation
    x_old((n+1)*(N+2)) = x((n+1)*(N+2));
end
% display(ii)

%% Write y
y = zeros(n,2);
for i = 1:n
    y(i,1) = x_old(1+(i-1)*(N+2));
    y(i,2) = x_old(i*(N+2));
end

% structure of y:
%     1    2
% 1 c1,m c1,p
% 2 c2,m c2,p
% ...
% n cn,m cn,p
end