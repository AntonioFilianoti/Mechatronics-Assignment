% direct_optimal_control_fmincon.m
% =========================================================================
% Direct optimal control (single shooting) using fmincon
% - Decision variables: discretized controls u1 and u2 (piecewise-linear)
% - Simulation: ode113 with griddedInterpolant for controls
% - Cost: integral running cost + terminal cost + soft obstacle penalty
% - No analytic gradients provided (fmincon approximates numerically)
% =========================================================================
clear; close all; clc;

%% ------------------------ Problem definition -----------------------------
m   = 1;       % massa [kg]
Cd  = 0.35;    % coefficiente di resistenza
t0  = 0; tf = 4;

% Stati: z = [x; y; theta; v]
z0 = [0; 0; pi/4; 0];     % stato iniziale
z_target = [10; 10; 0; 1];% stato desiderato a tf

% Cost weights
R = diag([1, 1]);   % controllo
Q = 0.1;             % penalizzazione sulla velocità (termine interno)
P = diag([200000,200000, 1, 1]); % penalizzazione finale

% Soft circular obstacle (penalized in cost)
xc = 3; yc = 5; r = 2;
alpha = 500000; sigma = 10;  % parametri soft penalty

%% -------------------- Discretizzazione dei comandi ----------------------
Nsegment = 60;                     % numero di nodi di controllo (più basso = più veloce)
Tu = linspace(t0, tf, Nsegment)';  % colonne

% decision vector: U = [u1_1 ... u1_N, u2_1 ... u2_N]'  (2N x 1)
nVar = 2 * Nsegment;

% initial guess (colonna)
u1_init = 5 * ones(Nsegment,1);  % prova con qualche spinta iniziale
u2_init = zeros(Nsegment,1);
U0 = [u1_init; u2_init];

% bounds on controls
u_max = 50; u_min = -50;
lb = [u_min*ones(Nsegment,1); u_min*ones(Nsegment,1)];
ub = [u_max*ones(Nsegment,1); u_max*ones(Nsegment,1)];

%% ----------------------- fmincon options --------------------------------
opts = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1e5, ...
    'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-4, ...
    'PlotFcn', {'optimplotfval','optimplotx'});

% If Parallel Toolbox available, user can enable:
opts.UseParallel = true;

%% ----------------------- Call to fmincon --------------------------------
problem.objective = @(U) cost_direct(U, Tu, z0, z_target, m, Cd, R, Q, P, xc, yc, r, alpha, sigma);
problem.x0 = U0;
problem.lb = lb;
problem.ub = ub;
problem.solver = 'fmincon';
problem.options = opts;

fprintf('Launching fmincon (direct shooting) with Nsegment = %d ...\n', Nsegment);
[Uopt, Jopt, exitflag, output] = fmincon(problem);

fprintf('Done: exitflag=%d, Jopt=%.6f\n', exitflag, Jopt);
disp(output);

%% ------------------ Reconstruct trajectory & plot -----------------------
u1_opt = Uopt(1:Nsegment);
u2_opt = Uopt(Nsegment+1:end);

% Create interpolants for final simulation
u1I = griddedInterpolant(Tu, u1_opt, 'linear', 'nearest');
u2I = griddedInterpolant(Tu, u2_opt, 'linear', 'nearest');

% simulate final
options = odeset('RelTol',1e-4,'AbsTol',1e-5);
[Tz, Z] = ode113(@(t,z) state_ode_for_direct(t,z,u1I,u2I,m,Cd), [t0 tf], z0, options);

% plot environment and trajectory
theta_plot = linspace(0,2*pi,300);
x_cons = xc + r*cos(theta_plot);
y_cons = yc + r*sin(theta_plot);

figure('Position',[200 200 900 650]);
plot(x_cons, y_cons, 'k--','LineWidth',1.5); hold on;
plot(Z(:,1), Z(:,2),'b-','LineWidth',2);
scatter(z0(1), z0(2), 80, 'g', 'filled');
scatter(z_target(1), z_target(2), 80, 'r', 'filled');
xlabel('x'); ylabel('y'); axis equal; grid on;
legend('Vincolo (soft)','Traiettoria ottima','Start','Goal');
title('Traiettoria ottenuta con Direct Shooting + fmincon');

% plot controls
figure('Position',[300 300 900 400]);
subplot(2,1,1); stairs(Tu, u1_opt, 'LineWidth',1.5); ylabel('u_1 = TPS'); grid on;
subplot(2,1,2); stairs(Tu, u2_opt, 'LineWidth',1.5); ylabel('u_2 = Steering Speed'); xlabel('t [s]'); grid on;

%% ===================== Local functions ==================================

function J = cost_direct(U, Tu, z0, z_target, m, Cd, R, Q, P, xc, yc, r, alpha, sigma)
% cost_direct - objective for direct shooting
%   INPUT:
%     U = [u1; u2] (2N x 1) decision vector (values at Tu nodes)
%     Tu: time nodes (N x 1)
%     z0, z_target: initial and desired final states
%     m, Cd: dynamics params
%     R, Q, P: cost weights
%     xc,yc,r,alpha,sigma: soft obstacle parameters
%   OUTPUT:
%     J: scalar cost (to minimize)
%
%  The cost is:
%    J = 0.5*(z(tf)-z_target)'*P*(z(tf)-z_target) + integral( 0.5*u'*R*u + Q*v^3 + soft_bound(x,y) dt )

    N = length(Tu);
    u1 = U(1:N);
    u2 = U(N+1:2*N);

    % create fast interpolants for controls (piecewise-linear)
    u1I = griddedInterpolant(Tu, u1, 'linear', 'nearest');
    u2I = griddedInterpolant(Tu, u2, 'linear', 'nearest');

    % simulate dynamics
    options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    [Tz, Z] = ode113(@(t,z) state_ode_for_direct(t,z,u1I,u2I,m,Cd), [Tu(1) Tu(end)], z0, options);

    ztf = Z(end,:)';

    % terminal cost
    J_terminal = 0.5 * (ztf - z_target)' * P * (ztf - z_target);

    % running cost: sample integrand at Tz and integrate (trap rule)
    % obtain u and v on Tz
    u1_on_Tz = u1I(Tz);
    u2_on_Tz = u2I(Tz);
    v_on_Tz  = Z(:,4);

    % control quadratic term
    uvecsq = 0.5*(R(1,1).*u1_on_Tz.^2 + R(2,2).*u2_on_Tz.^2); % assuming R diagonal
    % velocity penalty (example: Q * v^3)
    vterm = Q .* (v_on_Tz.^3);

    % soft boundary penalization
    x_on_Tz = Z(:,1);
    y_on_Tz = Z(:,2);
    expo = (r^2 - (x_on_Tz - xc).^2 - (y_on_Tz - yc).^2) ./ sigma;
    expo = min(expo, 700); % clamp to avoid overflow
    soft_term = alpha .* exp(expo);

    integrand = uvecsq + vterm + soft_term;

    % numerical integration (trapezoidal)
    J_running = trapz(Tz, integrand);

    % total cost
    J = J_terminal + J_running;
end


function dz = state_ode_for_direct(t, z, u1I, u2I, m, Cd)
% state equations used in direct method (uses interpolants for controls)
% z = [x; y; theta; v]
    u1 = u1I(t);
    u2 = u2I(t);

    xdot = z(4) * cos(z(3));
    ydot = z(4) * sin(z(3));
    thdot = u2;
    vdot = (1/m) * (u1 - Cd * z(4)^2);

    dz = [xdot; ydot; thdot; vdot];
end
