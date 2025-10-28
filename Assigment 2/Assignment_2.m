%% ========================================================================
%  Ottimizzazione del controllo di un veicolo planare con vincolo morbido
%  Descrizione:
%  Questo script risolve un problema di controllo ottimo per un veicolo
%  con dinamica non lineare, attrito aerodinamico e penalizzazione morbida
%  (soft constraint) su un ostacolo circolare nel piano (xc, yc, r).
%  ========================================================================

clear all; close all; clc;

%% ----------------------------- PARAMETRI BASE ----------------------------
m   = 1;        % Massa del veicolo [kg]
Cd  = 0.15;     % Coefficiente di resistenza aerodinamica
t0  = 0;        % Tempo iniziale [s]
tf  = 1.6;      % Tempo finale [s]

global nx nu N xc yc r h FigTraj

%% ----------------------- STATI E CONDIZIONI INIZIALI --------------------
% Stato: [x; y; theta; v]
% x, y   = coordinate posizione
% theta  = orientamento (rad)
% v      = velocità (m/s)

pos_i = [0; 0; 0];       % posizione e orientamento iniziali
pos_f = [1; 1; pi/3];    % posizione e orientamento finali

vel_i = 0;               % velocità iniziale
vel_f = 0;               % velocità finale

x_i = [pos_i; vel_i];    % stato iniziale
x_f = [pos_f; vel_f];    % stato finale desiderato

% Vettori di normalizzazione per scaling
xmax_vec = [pos_f(1:2); 2*pi/3; 1.5];  % valori massimi degli stati
umax_vec = [1; 1];                     % valori massimi dei controlli

%% -------------------------- PESI E VINCOLI COSTO -------------------------
w1 = 0.2;   % peso controllo 1
w2 = 0.3;   % peso controllo 2

% Peso finale sullo stato (errore finale rispetto a x_f)
P = diag(1 ./ xmax_vec.^2);
P = P .* diag([150 150 20 0]);  % penalizza fortemente errore su x, y, theta

% Peso sui controlli (energia spesa)
R = diag(1 ./ umax_vec.^2);
R = R .* diag([w1 w2]);

Q = 0.05;   % peso su penalizzazione della velocità

% Parametri del soft constraint (ostacolo circolare)
alpha = 2;       % peso del vincolo morbido
sigma = 0.01;    % "softness" → più piccolo = vincolo più rigido
xc = 0.65;       % centro ostacolo (x)
yc = 0.65;       % centro ostacolo (y)
r  = 0.15;       % raggio ostacolo

%% ----------------------------- DINAMICA -------------------------------
% Funzione di stato: x_dot = f(x,u)
% x = [x; y; theta; v], u = [F; omega]
% dove F è la forza di trazione e omega la variazione d’angolo (input angolare)
dx = @(x,u) [ ...
    x(4)*cos(x(3));                % x_dot = v*cos(theta)
    x(4)*sin(x(3));                % y_dot = v*sin(theta)
    u(2);                          % theta_dot = omega
    (1/m)*(u(1) - Cd*(x(4))^2) ];  % v_dot = (F - Cd*v^2)/m

% Jacobiano della dinamica rispetto allo stato (fx = ∂f/∂x)
fx = @(x,u) [ ...
     0,  0,  -x(4)*sin(x(3)),   cos(x(3));
     0,  0,   x(4)*cos(x(3)),   sin(x(3));
     0,  0,   0,                0;
     0,  0,   0,               -2*Cd*x(4)/m ];

% Jacobiano della dinamica rispetto al controllo (fu = ∂f/∂u)
fu = @(x,u) [ ...
     0,  0;
     0,  0;
     0,  1;
     1/m, 0];

%% -------------------------- FUNZIONE DI COSTO --------------------------
% Penalizzazione morbida vicino all’ostacolo (soft constraint)
Soft_cost_fun = @(x,u) alpha * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);

% Costo istantaneo (running cost)
L = @(x,u) 0.5*(u(1).^2*R(1,1) + u(2).^2*R(2,2)) ...  % costo controllo
          + Q*(x(4).^3) ...                           % penalizza velocità elevata
          + Soft_cost_fun(x,u);                       % penalità prossimità ostacolo

% Costo finale (terminal cost)
p = @(x) 0.5*(x - x_f)'*P*(x - x_f);

%% -------------------------- DERIVATE DEL COSTO -------------------------
% Gradienti di L rispetto a x e u (necessari per ottimizzazione)

a = @(x,u) -alpha * (2*(x(1)-xc)/sigma) * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);
b = @(x,u) -alpha * (2*(x(2)-yc)/sigma) * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);

Lx = @(x,u) [a(x,u), b(x,u), 0, 3*Q*x(4).^2];          % ∂L/∂x
Lu = @(x,u) [R(1,1)*u(1), R(2,2)*u(2)];                % ∂L/∂u

px = @(x) P*(x - x_f);                                 % ∂p/∂x

%% ------------------------- PARAMETRI DISCRETIZZAZIONE ------------------
nx = 4;                    % numero stati
nu = 2;                    % numero controlli
N  = 301;                  % numero di passi temporali
h  = tf / (N - 1);         % passo temporale (Δt)

%% -------------------- OPZIONI DELL'OTTIMIZZATORE ----------------------
options = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'Display', 'iter', ...
    'MaxIterations', 40, ...
    'OutputFcn', @plotTrajectoryCallback, ...
    'PlotFcn', [], ...
    'OptimalityTolerance', 1e-5, ...   % arresto se gradiente piccolo
    'StepTolerance', 1e-7, ...         % arresto se il passo è piccolo
    'FunctionTolerance', 1e-5);        % arresto se J varia poco

%% --------------------- STRUTTURA PARAMETRI PER FMINCON -----------------
param.N  = N;
param.nu = nu;
param.nx = nx;
param.dx = dx;
param.fx = fx;
param.fu = fu;
param.L  = L;
param.Lx = Lx;
param.Lu = Lu;
param.p  = p;
param.px = px;
param.h  = h;
param.x_i = x_i;
param.xc = xc;
param.yc = yc;
param.r  = r;

%% ---------------------- STIMA INIZIALE PER L’OTTIMIZZAZIONE ------------
% Costruzione di una traiettoria iniziale di riferimento
z0 = build_initial_guess(param, x_i, x_f, m, Cd);

%% --------------------- DEFINIZIONE PROBLEMA DI MINIMIZZAZIONE ----------
% Funzione obiettivo (costo totale e gradiente)
ObjFun = @(z) cost_and_grad(z, param);

% Vincoli non lineari (dinamica e ostacoli)
NLcon = @(z) con_and_grad(z, param);

% Vincoli lineari (qui assenti)
A = []; b = [];
Aeq = []; beq = [];
lb = []; ub = [];  % nessun limite sui controlli o stati

%% -------------------- RISOLUZIONE CON FMINCON --------------------------
[z, fval] = fmincon(ObjFun, z0, A, b, Aeq, beq, lb, ub, NLcon, options);

% z contiene la traiettoria ottima (stati + controlli)
% fval è il valore finale del costo

