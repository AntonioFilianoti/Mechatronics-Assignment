%% ========================================================================
%  Ottimizzazione MINIMUM TIME di un veicolo planare con vincolo morbido
%  ========================================================================
clear all; close all; clc;

%% ----------------------------- PARAMETRI BASE ----------------------------
m   = 1;        % Massa del veicolo [kg]
Cd  = 0.15;     % Coefficiente di resistenza aerodinamica
t0  = 0;        % Tempo iniziale [s]
% tf_guess ora è solo una stima iniziale per l'ottimizzatore
tf_guess  = 1.6;      

global nx nu N xc yc r
%% ----------------------- STATI E CONDIZIONI INIZIALI --------------------
% Stato: [x; y; theta; v]
pos_i = [0; 0; 0];       
pos_f = [1; 1; pi/3];    
vel_i = 0;               
vel_f = 0;               
x_i = [pos_i; vel_i];    
x_f = [pos_f; vel_f];    

% Limiti fisici (IMPORTANTE per Minimum Time)
v_max = 8;             % Velocità massima
u_max_force = 8.0;       % Forza massima
u_max_omega = 8.0;       % Velocità angolare massima

%% -------------------------- PESI E VINCOLI COSTO -------------------------
w1 = 0.2;   % peso controllo 1
w2 = 0.3;   % peso controllo 2

% Vettori di normalizzazione per scaling
xmax_vec = [pos_f(1:2); 2*pi/3; 1.5];  % valori massimi degli stati
umax_vec = [1; 1];                     % valori massimi dei controlli
% Peso finale sullo stato (errore finale rispetto a x_f)
P = diag(1 ./ xmax_vec.^2);
P = P .* diag([150 150 20 0]);  % penalizza fortemente errore su x, y, theta

% Peso sui controlli (energia spesa)
R = diag(1 ./ umax_vec.^2);
R = R .* diag([w1 w2]);

Q = 0.05;   % peso su penalizzazione della velocità

% Parametri del soft constraint (ostacolo circolare)
alpha = 12;       % peso del vincolo morbido
sigma = 0.01;    % "softness" → più piccolo = vincolo più rigido
xc = 0.65;       % centro ostacolo (x)
yc = 0.65;       % centro ostacolo (y)
r  = 0.15;       % raggio ostacolo


% Peso predominante sul tempo finale
Weight_Tf = 10.0; 

%% ----------------------------- DINAMICA -------------------------------
dx = @(x,u) [ ...
    x(4)*cos(x(3));                
    x(4)*sin(x(3));                
    u(2);                          
    (1/m)*(u(1) - Cd*(x(4))^2) ];  

fx = @(x,u) [ ...
     0,  0,  -x(4)*sin(x(3)),   cos(x(3));
     0,  0,   x(4)*cos(x(3)),   sin(x(3));
     0,  0,   0,                0;
     0,  0,   0,               -2*Cd*x(4)/m ];

fu = @(x,u) [ ...
     0,  0;
     0,  0;
     0,  1;
     1/m, 0];

%% -------------------------- FUNZIONE DI COSTO --------------------------
Soft_cost_fun = @(x,u) alpha * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);
% Gradienti costo soft
a = @(x,u) -alpha * (2*(x(1)-xc)/sigma) * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);
b = @(x,u) -alpha * (2*(x(2)-yc)/sigma) * exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2) / sigma);
Lx_soft = @(x,u) [a(x,u), b(x,u), 0, 0]; 

%% ------------------------- PARAMETRI DISCRETIZZAZIONE ------------------
nx = 4; nu = 2; N  = 101; % Nodi
% Il passo h non è calcolato qui, dipende dalla variabile decisionale tf

%% --------------------- STRUTTURA PARAMETRI -----------------
param.N  = N; param.nu = nu; param.nx = nx;
param.dx = dx; param.fx = fx; param.fu = fu;
param.Lx_soft = Lx_soft;
param.Soft_cost_fun = Soft_cost_fun;
param.R  = R;
param.Weight_Tf = Weight_Tf;
param.x_i = x_i; param.x_f = x_f;
param.xc = xc; param.yc = yc; param.r  = r;

%% ---------------------- STIMA INIZIALE E BOUNDS ------------------------
% Costruiamo z0. Nota: z0 deve avere lunghezza N*(nx+nu) + nx + 1 (per tf)
z0_traj = build_initial_guess(param, x_i, x_f, m, Cd, tf_guess);
z0 = [z0_traj; tf_guess]; % Aggiungiamo tf in fondo al vettore

% Definizione dei limiti (Bounds) - CRUCIALE per Minimum Time
% z = [x1, u1, x2, u2 ... xN, uN, xN+1, tf]
num_vars = numel(z0);
lb = -inf(num_vars, 1);
ub =  inf(num_vars, 1);

% Applichiamo i limiti ai controlli e velocità lungo tutta la traiettoria
for k = 1:N
    idx_u = (k-1)*(nx+nu) + nx + (1:nu);
    idx_v = (k-1)*(nx+nu) + 4; % Indice della velocità nello stato
    
    % Limiti controlli
    lb(idx_u) = [-u_max_force; -u_max_omega];
    ub(idx_u) = [ u_max_force;  u_max_omega];
    
    % Limite velocità (opzionale ma consigliato)
    lb(idx_v) = -v_max; 
    ub(idx_v) =  v_max; 
end
% Limite velocità stato finale
lb(end-nx) = -v_max; ub(end-nx) = v_max; % indice 4 dello stato finale (v) è (end-1-nx+4)

% Limiti per il tempo finale tf
idx_tf = num_vars; 
lb(idx_tf) = 0.5;   % Tempo minimo fisico plausibile
ub(idx_tf) = 10.0;  % Tempo massimo

%% -------------------- OPZIONI DELL'OTTIMIZZATORE ----------------------
options = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'Display', 'iter', ...
    'MaxIterations', 200, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-6, ...
    'CheckGradients', false); % Imposta true per debuggare i gradienti

%% -------------------- RISOLUZIONE CON FMINCON --------------------------
ObjFun = @(z) cost_and_grad_MinTime(z, param);
NLcon  = @(z) con_and_grad_MinTime(z, param);

[z_opt, fval] = fmincon(ObjFun, z0, [], [], [], [], lb, ub, NLcon, options);

tf_opt = z_opt(end);
fprintf('Tempo Minimo Ottenuto: %.4f s\n', tf_opt);

%% -------------------- PLOT RISULTATI --------------------------
plot_results(z_opt, param);

%% ========================================================================
%  FUNZIONI AUSILIARIE (Cost, Constraints, Guess, Plot)
%  ========================================================================

function [cost, grad] = cost_and_grad_MinTime(z, param)
    % Estrai tf (ultima variabile)
    tf = z(end);
    % Rimuovi tf per ottenere il vettore traiettoria classico
    z_traj = z(1:end-1);
    
    N  = param.N; nx = param.nx; nu = param.nu;
    h  = tf / (N - 1); % h è variabile ora!
    R  = param.R;
    Soft_cost = param.Soft_cost_fun;
    Lx_soft   = param.Lx_soft;
    W_tf = param.Weight_Tf;

    % Ricostruisci matrici
    x = zeros(nx, N+1); u = zeros(nu, N);
    for k = 1:N+1
        idx = (k-1)*(nx+nu) + (1:nx);
        x(:,k) = z_traj(idx);
    end
    for k = 1:N
        idx = (k-1)*(nx+nu) + nx + (1:nu);
        u(:,k) = z_traj(idx);
    end

    % --- CALCOLO COSTO ---
    running_cost = 0;
    soft_obs_cost = 0;
    
    % Costo = W_tf * tf + termine regolarizzazione controlli + ostacoli
    % Integrale L(x,u) dt approx sum(L * h)
    for k = 1:N
        ctrl_cost = 0.5*(u(:,k)' * R * u(:,k));
        obs_cost  = Soft_cost(x(:,k), u(:,k));
        running_cost = running_cost + ctrl_cost;
        soft_obs_cost = soft_obs_cost + obs_cost;
    end
    
    % Funzione Obiettivo Totale
    cost = W_tf * tf + h * (running_cost + soft_obs_cost);

    % --- CALCOLO GRADIENTE ---
    if nargout > 1
        grad = zeros(size(z));
        
        % Derivata parziale rispetto a tf (dJ/dtf)
        % dJ/dtf = W_tf + (running_cost + soft_obs_cost) * (dh/dtf)
        % con dh/dtf = 1/(N-1)
        dJ_dtf = W_tf + (running_cost + soft_obs_cost) / (N-1);
        
        % Derivate rispetto a stati e controlli
        for k = 1:N
            idx_x = (k-1)*(nx+nu) + (1:nx);
            idx_u = (k-1)*(nx+nu) + nx + (1:nu);
            
            % Gradiente costo ostacolo
            grad_obs = Lx_soft(x(:,k), u(:,k)); 
            
            % Contributo: h * dL/dx
            grad(idx_x) = grad(idx_x) + h * grad_obs'; 
            
            % Contributo: h * dL/du (qui solo R*u)
            grad(idx_u) = grad(idx_u) + h * (R * u(:,k));
        end
        
        % Assegna gradiente ultima componente
        grad(end) = dJ_dtf; 
    end
end

function [c, ceq, gradc, gradceq] = con_and_grad_MinTime(z, param)
    % Estrai tf e parametri
    tf = z(end);
    z_traj = z(1:end-1);
    
    N = param.N; nx = param.nx; nu = param.nu;
    h = tf / (N - 1);
    
    dx = param.dx; fx = param.fx; fu = param.fu;
    x_i = param.x_i; x_f = param.x_f;
    xc = param.xc; yc = param.yc; r = param.r;
    
    % Ricostruisci stati e controlli
    x = zeros(nx, N+1); u = zeros(nu, N);
    for k = 1:N+1
        x(:,k) = z_traj((k-1)*(nx+nu) + (1:nx));
    end
    for k = 1:N
        u(:,k) = z_traj((k-1)*(nx+nu) + nx + (1:nu));
    end
    
    % --- VINCOLI UGUAGLIANZA (Dinamica + Boundary) ---
    % Dinamica: x(k+1) - x(k) - h*f(x,u) = 0
    % Init: x(1) - x_i = 0
    % Final: x(N+1) - x_f = 0  (Vincolo hard sullo stato finale)
    
    num_eq = nx*N + nx + nx; % Dinamica + Init + Final
    ceq = zeros(num_eq, 1);
    
    % Inizializza Jacobiano sparso
    nz = numel(z);
    gradceq = sparse(nz, num_eq);
    I = eye(nx);
    
    % 1. Dinamica
    for k = 1:N
        idx_eq = (k-1)*nx + (1:nx);
        fk = dx(x(:,k), u(:,k));
        
        % Vincolo
        ceq(idx_eq) = x(:,k+1) - x(:,k) - h*fk;
        
        % Derivate rispetto a z_traj
        Ak = -I - h*fx(x(:,k), u(:,k)); % d/dx_k
        Bk = -h*fu(x(:,k), u(:,k));     % d/du_k
        Ck =  I;                        % d/dx_{k+1}
        
        zxk   = (k-1)*(nx+nu) + (1:nx);
        zuk   = (k-1)*(nx+nu) + nx + (1:nu);
        zxkp1 = (k)*(nx+nu) + (1:nx);
        
        gradceq(zxk,   idx_eq) = Ak';
        gradceq(zuk,   idx_eq) = Bk';
        gradceq(zxkp1, idx_eq) = Ck';
        
        % Derivata rispetto a tf (Nuovo!)
        % d/dtf (ceq) = - f(x,u) * dh/dtf = - f(x,u) * 1/(N-1)
        gradceq(end, idx_eq) = (-fk / (N-1))'; 
    end
    
    % 2. Condizione Iniziale
    idx_eq_i = nx*N + (1:nx);
    ceq(idx_eq_i) = x(:,1) - x_i;
    gradceq(1:nx, idx_eq_i) = I;
    
    % 3. Condizione Finale (Hard Constraint)
    idx_eq_f = nx*N + nx + (1:nx);
    ceq(idx_eq_f) = x(:,end) - x_f;
    idx_x_end = N*(nx+nu) + (1:nx);
    gradceq(idx_x_end, idx_eq_f) = I;
    
    % --- VINCOLI DISUGUAGLIANZA (Ostacolo) ---
    c = zeros(N+1, 1);
    gradc = sparse(nz, N+1);
    
    for k = 1:N+1
        dist_sq = (x(1,k)-xc)^2 + (x(2,k)-yc)^2;
        c(k) = r^2 - dist_sq; % c <= 0 se fuori dall'ostacolo
        
        % Gradiente
        if c(k) > -1e-2 % Ottimizzazione: calcola gradiente solo se vicino
            gx = -2*(x(1,k)-xc);
            gy = -2*(x(2,k)-yc);
            idx_x = (k-1)*(nx+nu) + (1:nx); % attenzione a ultimo stato
            if k == N+1, idx_x = N*(nx+nu) + (1:nx); end
            
            gradc(idx_x(1), k) = gx;
            gradc(idx_x(2), k) = gy;
        end
    end
end

function z0 = build_initial_guess(param, x_i, x_f, m, Cd, tf)
    % Usa interpolazione semplice
    N = param.N; nx = param.nx; nu = param.nu; h = tf/(N-1);
    x = zeros(nx, N+1); u = zeros(nu, N);
    
    % Lineare da start a end
    for i=1:4
        x(i,:) = linspace(x_i(i), x_f(i), N+1);
    end
    
    % Guess controlli (molto approssimato)
    dx_diff = diff(x(1,:)); dy_diff = diff(x(2,:));
    vel = sqrt(dx_diff.^2 + dy_diff.^2)/h;
    x(4,1:end-1) = vel; x(4,end) = 0;
    
    % Impacchettamento
    z0 = zeros(N * (nx + nu) + nx, 1);
    for k = 1:N
        base = (k - 1) * (nx + nu);
        z0(base + (1:nx))      = x(:, k);
        z0(base + nx + (1:nu)) = u(:, k);
    end
    z0(N * (nx + nu) + (1:nx)) = x(:, end);
end

function plot_results(z, param)
    tf = z(end);
    z_traj = z(1:end-1);
    N = param.N; nx = param.nx; nu = param.nu; h = tf/(N-1);
    
    x = zeros(nx, N+1); u = zeros(nu, N);
    for k = 1:N+1
        x(:,k) = z_traj((k-1)*(nx+nu) + (1:nx));
    end
    for k = 1:N
        u(:,k) = z_traj((k-1)*(nx+nu) + nx + (1:nu));
    end
    t = 0:h:tf;
    
    figure('Color','w');
    subplot(2,2,[1 3]); hold on; grid on; axis equal;
    viscircles([param.xc param.yc], param.r, 'Color', 'r');
    plot(x(1,:), x(2,:), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    quiver(x(1,1:10:end), x(2,1:10:end), cos(x(3,1:10:end)), sin(x(3,1:10:end)), 0.1, 'k');
    xlabel('X [m]'); ylabel('Y [m]'); title(['Traiettoria Ottima (T_f = ' num2str(tf, '%.3f') ' s)']);
    
    subplot(2,2,2); 
    stairs(t, u(1,:), 'LineWidth', 1.5); grid on;
    ylabel('Forza [N]'); title('Controllo u1 (Trazione)');
    
    subplot(2,2,4); 
    stairs(t, u(2,:), 'LineWidth', 1.5); grid on;
    xlabel('Tempo [s]'); ylabel('Omega [rad/s]'); title('Controllo u2 (Sterzo)');
end