% Calculate the required pressure to for the given membrane
function output = solvePressure(membrane_paras, Q0, temp, c0, WR_target, moduleParas, calParas)

rp = membrane_paras(1);
xe = membrane_paras(2);
A = (rp*10^9)^2*10^2*3.6*10^-9/xe/(8*8.9e-4);
S = moduleParas.N_leaf*moduleParas.W_module*moduleParas.L_module; % Total membrane area, m2

% Define 3 initial trial values of required pressure
c0_Mg = [0;c0(2);2*c0(2)]; % mM
delta_pi_f = 8.314*temp*sum(c0_Mg); % Pa
rMg = [0.9 0.5 0.2];
deltaP = delta_pi_f.*rMg.*(1-WR_target).^(-rMg);

solveWR = @(pressure) compute_WR_difference(pressure, membrane_paras, c0, temp, Q0, WR_target, moduleParas, calParas);

% Try P_0 in turn
P0_candidates = deltaP;
nCand = numel(P0_candidates);
failure_msgs = cell(nCand, 1);
found = false;
for k = 1:nCand
    P_try = P0_candidates(k);
    y = solveWR(P_try);
    if isfield(y, 'exitFlag') && y.exitFlag
        P_0 = P_try;
        WR_0 = y.results.WR;
        found = true;
        break;
    else
        if isfield(y, 'errorMessage') && isfield(y.errorMessage, 'message')
            failure_msgs{k} = y.errorMessage.message;
        else
            failure_msgs{k} = 'No errorMessage returned';
        end
    end
end

if ~found
    msg = sprintf('All initial pressure candidates cannot be used to calculate WR_0: \n');
    for k = 1:nCand
        msg = [ msg, ...
            sprintf('Candidate %d: P0 = %.2f bar — %s\n', ...
            k, P0_candidates(k)/1e5, failure_msgs{k}) ];
    end
    error(msg);
end

% Find root interval
[Plow, Phigh, nEval] = bracketPressure(solveWR, P_0, WR_0, WR_target);

% Use fzero to solve the pressure and record other results simultaneously
finalOutput = [];

    function y = WRdiff(pressure)
        out = compute_WR_difference(pressure, membrane_paras, c0, temp, Q0, WR_target, moduleParas, calParas);
        finalOutput = out;
        y = out.WR_difference;
    end

options = optimset('Display','off','TolX',1e-4);
[P_solution, fval, exitFlag, logs] = fzero(@WRdiff, [Plow, Phigh], options);

% Output
output.solvedP = P_solution;
output.fval = fval;
output.exitFlag = exitFlag;
output.fzeroLog = logs;
output.results = finalOutput.results;
end


%------------------------------------------------------------------

function output = compute_WR_difference(pressure, membrane_paras, c0, temp, Q0, WR_target, moduleParas, calParas)

output.exitFlag = false;

try
    % Define parameters
    rp = membrane_paras(1);
    xe = membrane_paras(2);
    n = length(c0);
    epsilon_b = 78.4;

    % Calibrating ion diffusivity at given temperature
    vis_u = exp(-10.5+530/(temp-146)); % μ, Pa·s

    % Module parameters
    N_leaf = moduleParas.N_leaf; % Number of leaves in each membrane module
    W_module = moduleParas.W_module; % Effective leaf width, m
    L_module = moduleParas.L_module; % Effective leaf length, m

    % Calculation options
    eta = calParas.eta; % Under-relaxation factor
    tolerance = calParas.tolerance; % Maximum normalized residuals
    N = calParas.N; % Number of discretization points in coupon scale
    dL = calParas.dL; % Step length in a module, m

    %------------------------------------------------------------------

    m = L_module/dL; % number of steps

    jvlist = zeros(m,1); % Local water flux, m/s

    Qpmlist = zeros(m+1,1); % Cumulative permeate flux rate, m3/s
    Mpmlist = zeros(m+1,n); % Cumulative permeate mass (mole) flux rate, mol/s
    wrlist = zeros(m,1); % Accumulative water recovery of step m

    cflist = zeros(m+1,n); % Local feed concentration, mol/m3
    cflist_cp = zeros(m,n); % Local feed concentration w/ cp, mol/m3
    cplist = zeros(m,n); % Local permeate concentration, mol/m3
    cpmlist = zeros(m,n); % Cumulative permeate concentration, mol/m3
    rlilist = zeros(m,1); % Li recovery
    purity = zeros(m,1); % Li purity
    kcs = zeros(m,n); % Mass transfer coefficient of each step, m/s

    cflist(1,:) = c0';

    operation_paras = {NaN, temp, NaN}'; % jv, c0 = NaN
    other_paras = {NaN, eta, tolerance, N}'; % kc = NaN

    eta_jv = 0.5;

    for i = 1:m
        if i == 1
            jv_guess = 1e-5; % m/s
            jv_old = 0;
            jv_limit = 0.001e-5;
        else
            jv_guess = 0.9*jv_guess;
            jv_limit = 0.001*jvlist(i-1);
        end

        j = 1;
        while abs(jv_guess-jv_old)>jv_limit
            if j>200 && eta_jv == 0.5
                eta_jv = 0.01;
                j = 1;
            end

            if j>300 && eta_jv == 0.01
                eta_jv = 0.005;
                j = 1;
            end

            if j>300 && eta_jv == 0.005
                error('jv_guess failed to converge (j>300). eta_jv = %g', eta_jv);
            end

            operation_paras{1} = jv_guess;

            % Calculating K_·c,i
            if i == 1
                WR_now = 0;
            else
                WR_now = wrlist(i-1);
            end
            kc_conv = solve_kc(N_leaf,W_module,Q0,WR_now,temp);
            kc = solve_kc_alt(kc_conv, jv_guess);
            other_paras{1} = kc;

            % Input local feed concentration into model
            operation_paras{3} = cflist(i,:)';

            % Call up DSPM-DE model
            fallbackEtas = [0.01, 0.005, 0.001];
            eta0 = other_paras{2};
            eta_candidates = [ eta0, fallbackEtas(fallbackEtas < eta0) ];

            success = false;
            lastME  = [];

            for eta_try = eta_candidates
                other_paras{2} = eta_try;
                try
                    y = model(membrane_paras, operation_paras, other_paras);
                    success = true;
                    break;
                catch ME
                    if strcmp(ME.message, 'Failed to converge, reduce eta.')
                        continue;
                    else
                        lastME = ME;
                        break;
                    end
                end
            end

            if ~success
                if isempty(lastME)
                    error('The model failed to converge for all values ​​of η [%s].', ...
                        num2str(eta_candidates));
                else
                    error('The model encountered an error: %s', lastME.message);
                end
            end

            cm = y(:,1); % mol/m3
            cp = y(:,2); % mol/m3

            pi_m = 8.314*temp*sum(cm); % Pa
            pi_p = 8.314*temp*sum(cp); % Pa
            delta_pi = pi_m-pi_p; % Pa

            jv_old = jv_guess; % m/s
            jv_guess = rp^2/(8*vis_u*xe)*(pressure-delta_pi); % m/s

            delta_jv = min(abs(jv_old/(jv_guess-jv_old)),1)*(jv_guess-jv_old)/jv_old;
            jv_guess = jv_old*(1+eta_jv*delta_jv);

            if jv_guess<0
                error('jv_guess<0');
            end

            % display(j)
            j = j+1;
        end

        jvlist(i) = jv_guess;
        Qpmlist(i+1) = Qpmlist(i)+jv_guess*N_leaf*W_module*dL; % m3/s
        if Qpmlist(i+1)>=Q0
            error('Q0 is too small');
        end
        wrlist(i) = Qpmlist(i+1)/Q0;

        cplist(i,:) = cp'; % mol/m3 = mM
        cflist_cp(i,:) = cm';
        Mpmlist(i+1,:) = Mpmlist(i,:)+jv_guess*N_leaf*W_module*dL*cplist(i,:); % mol/s
        rlilist(i) = Mpmlist(i+1,1)/(c0(1)*Q0);
        cpmlist(i,:) = Mpmlist(i+1,:)./Qpmlist(i+1); % mol/m3 = mM
        purity(i) = cpmlist(i,1)*6.94/(cpmlist(i,1)*6.94+cpmlist(i,2)*24.31);
        cflist(i+1,:) = ((Q0-Qpmlist(i)).*cflist(i,:)-jv_guess*N_leaf*W_module*dL.*cplist(i,:))/(Q0-Qpmlist(i+1));
        kcs(i,:) = kc';
    end

    summary.WR = wrlist(m);
    summary.LI_REC = rlilist(m);
    summary.LI_PUR = purity(m);
    summary.CC = cflist(m+1,:);
    summary.CP = cpmlist(m,:);

    summary.JV = jvlist;
    summary.CF_LOC = cflist;
    summary.CM_LOC = cflist_cp;
    summary.CP_LOC = cplist;
    summary.CP_CUMLIST = cpmlist;
    summary.LI_RECLIST = rlilist;
    summary.LI_PURLIST = purity;
    summary.QP_CUMLIST = Qpmlist(2:end);
    summary.MP_CUMLIST = Mpmlist(2:end,:);
    summary.WR_LIST = wrlist;
    summary.KC_LIST = kcs;

    output.WR_difference = wrlist(m)-WR_target;
    output.results = summary;
    output.exitFlag = true;
catch ME
    output.errorMessage = ME;
end
end


%------------------------------------------------------------------

function [Plow, Phigh, nEval] = bracketPressure(solveWR, P0, WR0, WR_target)
maxExpand = 6;
factorInc = 2;
maxBisect = 8;

nEval = 0;

if WR0 < WR_target
    % Probe toward the high-pressure side
    Psucc  = P0;
    WRsucc = WR0;
    step   = 0.5 * P0;    % Initial step size
    hitFail     = false;
    Pfail       = NaN;

    % (1) Geometric extrapolation
    for k = 1:maxExpand
        Ptrial = Psucc + step;
        output = solveWR(Ptrial);
        ok = output.exitFlag;
        if ok
            WRtrial = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRtrial > WR_target
                % Bracket found
                Plow  = Psucc;
                Phigh = Ptrial;
                return;
            else
                % Update successful point and increase step size
                Psucc  = Ptrial;
                WRsucc = WRtrial;
                step   = step * factorInc;
            end
        else
            % Record failure wall
            Pfail   = Ptrial;
            hitFail = true;
            break;
        end
    end

    if ~hitFail
        error('After %d extrapolations toward the high-pressure side, still no bracket or failure. Please increase maxExpand.', maxExpand);
    end

    % (2) Bisection fallback: within [Psucc, Pfail]
    a   = Psucc;
    b   = Pfail;

    for j = 1:maxBisect
        Pmid = 0.5 * (a + b);
        output = solveWR(Pmid);
        ok = output.exitFlag;
        if ok
            WRmid = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRmid > WR_target
                Plow  = a;
                Phigh = Pmid;
                return;
            else
                a = Pmid;
            end
        else
            b = Pmid;
        end
    end

    error('After %d bisections still no bracket found. Please increase maxBisect.', maxBisect);

else
    % Direct bisection toward the low-pressure side [0, P0]
    a   = 0;
    b   = P0;
    WRb = WR0;

    for j = 1:maxBisect
        Pmid = 0.5 * (a + b);
        output = solveWR(Pmid);
        ok = output.exitFlag;
        if ok
            WRmid = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRmid < WR_target
                % Bracket found
                Plow  = Pmid;
                Phigh = b;
                return;
            else
                b   = Pmid;
                WRb = WRmid;
            end
        else
            a = Pmid;
        end
    end

    error('After %d bisections toward the low-pressure side within [0, P0], still no bracket found. Please increase maxBisect.', maxBisect);
end
end