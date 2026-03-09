function output = solveFeed(membrane_paras, Q0, temp, pressure, WR_target, moduleParas, calParas)

rp = membrane_paras(1);
xe = membrane_paras(2);
A = (rp*10^9)^2*10^2*3.6*10^-9/xe/(8*8.9e-4);
S = moduleParas.N_leaf*moduleParas.W_module*moduleParas.L_module; % Total membrane area, m2
MLR = calParas.MLR; % The molar ratio of Mg to Li
NLR = calParas.NLR; % The molar ratio of Na to Li

% Guess initial c0 (Δπ = cf,MgCl2)
deltaP = Q0*3600*1000/(S*A)*log(1/(1-0.5)); % bar
delta_pi = pressure-deltaP*10^5; % Pa
c0_Mg_guess = delta_pi/(3*8.314*temp); % mM
c0_Li_guess = c0_Mg_guess/MLR;
c0_Na_guess = c0_Li_guess*NLR;
c0_guess = [c0_Li_guess; c0_Mg_guess; c0_Na_guess; c0_Li_guess+c0_Na_guess+2*c0_Mg_guess];
% c0_guess = [57.78	335.06	577.81	1305.70]';

solveWR = @(c0) compute_WR_difference(pressure, membrane_paras, c0, temp, Q0, WR_target, moduleParas, calParas);

y = solveWR(c0_guess);

if y.exitFlag
    WR_0 = y.results.WR;
else
    error('c0_guess error, msg: %s\n', y.errorMessage.message)
end

% Find root interval
[c0_low, c0_high, nEval] = bracketFeed(solveWR, c0_guess, WR_0, WR_target);

% Use fzero to solve the feed concentration and record other results simultaneously
finalOutput = [];

    function y = WRdiff(c0_Mg)
        c0_Li = c0_Mg/MLR;
        c0_Na = c0_Li*NLR;
        c0 = [c0_Li; c0_Mg; c0_Na; c0_Mg*2+c0_Li+c0_Na];
        out = compute_WR_difference(pressure, membrane_paras, c0, temp, Q0, WR_target, moduleParas, calParas);
        finalOutput = out;
        y = out.WR_difference;
    end

options = optimset('Display','off','TolX',1e-4);
[c0_Mg_solution, fval, exitFlag, logs] = fzero(@WRdiff, [c0_low(2), c0_high(2)], options);
c0_Li_solution = c0_Mg_solution/MLR;
c0_Na_solution = c0_Li_solution*NLR;
c0_solution = [c0_Li_solution; c0_Mg_solution; c0_Na_solution; c0_Mg_solution*2+c0_Li_solution+c0_Na_solution];

% Output
output.solvedFeed = c0_solution;
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

    % Calibrating ion diffusivity at given temperature
    vis_u = exp(-10.5+530/(temp-146));

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

function [c0_low, c0_high, nEval] = bracketFeed(solveWR, c0_guess, WR0, WR_target)
maxExpand = 6;
factorInc = 2;
maxBisect = 8;

nEval = 0;

if WR0 > WR_target
    c0_succ  = c0_guess;
    WRsucc = WR0;
    step   = 0.5 .* c0_guess;
    hitFail     = false;
    c0_fail       = NaN;

    % (1) Geometric extrapolation
    for k = 1:maxExpand
        c0_trial = c0_succ + step;
        output = solveWR(c0_trial);
        ok = output.exitFlag;
        if ok
            WRtrial = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRtrial < WR_target
                c0_low  = c0_succ;
                c0_high = c0_trial;
                return;
            else
                c0_succ  = c0_trial;
                WRsucc = WRtrial;
                step   = step * factorInc;
            end
        else
            c0_fail   = c0_trial;
            hitFail = true;
            break;
        end
    end

    if ~hitFail
        error('After %d extrapolations towards the higher concentration side, still no bracket or failure. Please increase maxExpand.', maxExpand);
    end

    % (2) Bisection fallback: within [Psucc, Pfail]
    a   = c0_succ;
    b   = c0_fail;

    for j = 1:maxBisect
        c0_mid = 0.5 * (a + b);
        output = solveWR(c0_mid);
        ok = output.exitFlag;
        if ok
            WRmid = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRmid < WR_target
                c0_low  = a;
                c0_high = c0_mid;
                return;
            else
                a = c0_mid;
            end
        else
            b = c0_mid;
        end
    end

    error('After %d bisections still no bracket found. Please increase maxBisect.', maxBisect);

else
    a   = 0;
    b   = c0_guess;
    WRb = WR0;

    for j = 1:maxBisect
        c0_mid = 0.5 * (a + b);
        output = solveWR(c0_mid);
        ok = output.exitFlag;
        if ok
            WRmid = output.results.WR;
        end
        nEval = nEval + 1;

        if ok
            if WRmid > WR_target
                c0_low  = c0_mid;
                c0_high = b;
                return;
            else
                b   = c0_mid;
                WRb = WRmid;
            end
        else
            a = c0_mid;
        end
    end

    error('After %d bisections on the lower concentration side within [0, c0_guess], the root was still not bracketed. Please increase maxBisect.', maxBisect);
end
end