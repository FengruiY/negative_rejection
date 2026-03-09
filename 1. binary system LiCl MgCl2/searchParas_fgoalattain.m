% Code for determining membrane parameters to achieve certain Mg2+ and Li+ rejection
% The results is used to calculate the SECLi reduction for NF step

clc; clear; tic
% Input
epsilon_p = 42.2; % Dielectric constant (i.e. relative permittivity) inside pore
epsilon_b = 78.4; % Relative permittivity for water

% Operation parameters
Q0 = 10/3600; % (m3/h)/3600, m3/s
temp = 298; % K
n = 3; % Specify the number of ions in the system
c0 = zeros(n,1);
% Input ion concentrations below
c0(1) = 4.611; % Li, mol/m3
c0(2) = 26.738; % Mg, mol/m3
c0(3) = c0(1)+2*c0(2); % Cl, mol/m3

% Module parameters
% 8040: D=200 mm, L=1 m
moduleParas.N_leaf = 20; % Number of leaves in each membrane module
moduleParas.W_module = 1*2; % Effective leaf width, m
moduleParas.L_module = 1*7; % Effective leaf length, m

% Calculation options
calParas.eta = 0.05; % Under-relaxation factor
calParas.tolerance = 1e-4; % Maximum normalized residuals
calParas.N = 50; % Number of discretization points in coupon scale
calParas.dL = 0.1; % Step length in a module, m

% ----------------------------------- Input fixed parameters above ---------------------------------------
% Set lower and upper bounds
lb = [0.351e-9, 0.5]; % rp (m), A (LMH/bar)
ub = [0.600e-9, 100];
x0 = [0.4e-9, 10];

% Load targets
T = readtable('caseList.xlsx');
nCase = height(T);
RMgList = T.RMg;
RLiList = T.RLi;
rhoMgList = T.rhoMg;
rhoLiList = T.rhoLi;
WRList = T.WR;
cxList = T.cx;

% Pre-define result storage
results = cell(nCase,1);

parfor k = 1:nCase
    % Generate logs for calculation progress
    logName = sprintf("case_%02d_log.txt", k);
    fid     = fopen(logName, "w");
    cleanupObj = onCleanup(@() fclose(fid));

    % fgoalattain options
    optFG = optimoptions('fgoalattain', ...
        'Display', 'none', ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', 1e-3, ...
        'StepTolerance', 1e-4, ...
        'FunctionTolerance', 1e-4, ...
        'ConstraintTolerance', 0.02, ...
        'OptimalityTolerance', 1e-4, ...
        'MaxIterations', 20, ...
        'UseParallel', false, ...
        'FunValCheck', 'on', ...
        'OutputFcn', @(x,ov,state) iterLogFcn(x,ov,state,fid) ...
        );
    
    % Read a target
    goal_Mg = RMgList(k);
    goal_Li = RLiList(k);
    rhoMg   = rhoMgList(k);
    rhoLi   = rhoLiList(k);

    goal   = [  goal_Mg,    goal_Li,   -goal_Mg,    -goal_Li ];
    weight = [   rhoMg,      rhoLi,     rhoMg,       rhoLi   ];

    WR_target = WRList(k);
    cx = cxList(k);

    try
        res = runCase(WR_target, cx, Q0, temp, c0, epsilon_p, moduleParas, calParas, lb, ub, x0, goal, weight, optFG, fid);
    catch ME
        res = struct();
        res.WR    = WR_target;
        res.cx    = cx;
        res.error = ME.message;
        fprintf('Error at case %d (WR=%g, cx=%g): %s\n', k, WR_target, cx, ME.message);
    end
    results{k} = res;
    disp(k)
end

% Save results
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = ['searchParameters_', timestamp, '.mat'];
save(filename);
fprintf('Results saved to %s\n', filename);

toc

%------------------------------------------------------------------

function f = runCase(WR_target, cx, Q0, temp, c0, epsilon_p, moduleParas, calParas, lb, ub, x0, goal, weight, optFG, fid)

lastOut = [];
nCall = 0;

% Normalize x0
z0 = (x0 - lb) ./ (ub - lb);

% Target function
vecFun = @(z) objVec(z);
    function y = objVec(z)
        rp = lb(1) + z(1) * (ub(1) - lb(1));
        A  = lb(2) + z(2) * (ub(2) - lb(2));
        xe = (rp*10^9)^2*10^2*3.6/(A*8*8.9e-4)*10^-9;

        key   = [rp, A];

        nCall = nCall + 1;

        % Check if this parameter pair has been tried
        [hit, yCache, outCache] = cacheManager('get', key);
        if hit
            y       = yCache;
            lastOut = outCache;
            fprintf(fid,"#%d | rp = %.5g  A = %.5g  (cached)\n", nCall, rp*10^9, A);
            return
        end

        fprintf(fid,"#%d | rp = %.5g  A = %.5g\n", nCall, rp*10^9, A);
        
        membrane_paras = [rp, xe, cx, epsilon_p]';
        outLocal = solvePressure(membrane_paras, Q0, temp, c0, WR_target, moduleParas, calParas);

        lastOut = outLocal;

        R = 1 - outLocal.results.CP ./ c0';
        R_Mg = R(2);
        R_Li = R(1);

        y = [R_Mg; R_Li; -R_Mg; -R_Li];

        % Update current result to cache
        cacheManager('add', key, y, outLocal);
    end

z_lb = [0, 0];
z_ub = [1, 1];

% Call up fgoalattain
[zSol, ~, attainfactor, flag, outFG] = fgoalattain( ...
    vecFun , z0 , goal , weight , ...
    [] , [] , [] , [] , ...       % A, b, Aeq, beq
    z_lb , z_ub , [] , ...            % lb, ub, nonlcon
    optFG );

xSol = lb + zSol .* (ub - lb);

% Output results
f.WR       = WR_target;
f.cx       = cx;
f.rp       = xSol(1);
f.A        = xSol(2);
f.gamma    = attainfactor;
f.exitFlag = flag;
f.logs     = outFG;
f.out_all  = lastOut; % save all other outputs of solvePressure
end


%------------------------------------------------------------------
% Function to write logs
function stop = iterLogFcn(~, ov, state, fid)
if strcmp(state, "iter")
    fvals = ov.fval;
    fvalStr = sprintf(' %.3g', fvals);

    fprintf(fid, ...
        "Iter=%3d  Fcount=%4d  γ=%8.3g  fval:%s  MaxConstr=%8.3g  step=%8.3g  1stOpt=%8.3g\n", ...
        ov.iteration, ov.funccount, ov.attainfactor, fvalStr, ...
        ov.constrviolation, ov.stepsize, ov.firstorderopt);

    cacheManager('reset');
elseif strcmp(state,"init")
    cacheManager('reset');
end
stop = false;
end


%------------------------------------------------------------------
% Cache recent trial results
function varargout = cacheManager(cmd, key, y, outLocal)
persistent keys ys outs
if isempty(keys)
    keys = zeros(0,2);  ys = {};  outs = {};
end
tol = 1e-14;

switch cmd
    case 'get'
        idx = find(all(abs(keys - key) < tol,2), 1);
        if isempty(idx)
            varargout = {false, [], []};
        else
            varargout = {true, ys{idx}, outs{idx}};
        end

    case 'add'
        keys = [keys; key];
        ys{end+1}   = y;
        outs{end+1} = outLocal;
        varargout = {};

    case 'reset'
        keys = zeros(0,2);
        ys   = {};
        outs = {};
        varargout = {};

    otherwise
        error("cacheManager: unknown cmd '%s'", cmd);
end
end