%% ========================================================================
%  Optimal Control with Adjoint Equations (Nonlinear System)
%  ------------------------------------------------------------------------
%  Questo script risolve un problema di controllo ottimo per un sistema 
%  dinamico non lineare (con attrito quadratico e vincolo morbido di traiettoria)
%  utilizzando un approccio iterativo basato su equazioni aggiunte (costate).
%  ------------------------------------------------------------------------
%% ========================================================================

clear; close all; clc;

%% ------------------------ Parametri del modello -------------------------
global m Cd mu_r epsilon
m   = 1;          % massa del corpo [kg]
Cd  = 0.15;       % coefficiente di resistenza aerodinamica (attrito quadratico)
epsilon = 1;
mu_r = 1;
t0  = 0;          % tempo iniziale [s]
tf  = 4;        % tempo finale [s]

% Stati del sistema: [x; y; theta; v]
pos_i = [0; 0; 0];         % posizione e orientamento iniziali (x, y, theta)
pos_f = [2; 2; pi/3];      % posizione e orientamento finali desiderati

vel_i = 0;                 % velocità iniziale
vel_f = 0;                 % velocità finale desiderata

z_i = [pos_i; vel_i];      % stato iniziale completo
z_f = [pos_f; vel_f];      % stato finale completo

% Vettori massimi per la normalizzazione degli stati e dei controlli
xmax_vec = [pos_f(1:2); 2*pi/3; 1.5];   % massimi per gli stati
umax_vec = [1, 1]';                     % massimi per i controlli

%% ----------------------- Pesi e vincoli del costo -----------------------
w1 = 0.2;                 % peso sul controllo 1
w2 = 0.2;                 % peso sul controllo 2

% Peso sullo stato finale (matrice P)
P = diag(1./xmax_vec.^2);
P = P .* diag([155 155 20 0]);  % maggior peso su x, y, theta; minore su v

% Peso sui controlli (matrice R)
R = diag(1./umax_vec.^2);
R = R .* diag([w1 w2]);

Q = 0.05;                 % peso sulla penalizzazione della velocità

% Parametri del vincolo morbido (soft constraint)
alpha = 2;                % peso della penalizzazione
sigma = 0.01;             % parametro di “morbidezza” (quanto il vincolo è sfumato)
xc = 0.65; yc = 0.65;     % coordinate del centro del vincolo
r = 0.2;                  % raggio della regione di vincolo

%% -------------------------- Plot del vincolo ----------------------------
theta = linspace(0, 2*pi, 400);        % angolo per disegnare la circonferenza
x_cons = xc + r .* cos(theta);         % coordinate x del vincolo
y_cons = yc + r .* sin(theta);         % coordinate y del vincolo

figure('Position', [100 100 900 700]);
plot(x_cons, y_cons, 'k--', 'LineWidth', 1.5); hold on;     % vincolo morbido
scatter(pos_i(1), pos_i(2), 70, 'g', 'filled');             % punto iniziale
scatter(pos_f(1), pos_f(2), 70, 'r', 'filled');             % punto finale
xlabel('x'); ylabel('y');
axis equal;
grid on; grid minor;
legend('Vincolo morbido', 'Start', 'Goal');
title('Problema di controllo ottimo con vincolo morbido');

%% --------------------------- Discretizzazione ---------------------------
Nsegment = 2000;                       % numero di intervalli temporali
Tu = linspace(t0, tf, Nsegment);       % vettore del tempo per i controlli
options = odeset('RelTol',1e-4,'AbsTol',1e-4);              % tolleranze ODE per costate
options_state = odeset('RelTol',1e-12,'AbsTol',1e-14);      % tolleranze ODE per stati

%% ----------------------- Parametri dell’iterazione ----------------------
Nmax = 5*10000;                         % numero massimo di iterazioni
step = 5e-3;                            % passo di aggiornamento del controllo
eps  = 0.1;                            % soglia di convergenza
u = [ones(1, Nsegment); zeros(1, Nsegment)];  % controllo iniziale (u1=1, u2=0)

%% ------------------------- Procedura Iterativa --------------------------
for ii = 1:Nmax

    % Integrazione delle equazioni di stato in avanti nel tempo
    [Tz, Z] = ode45(@(t,z) stateEq_A3(t,z,u,Tu,m ,Cd, mu_r, epsilon), [t0 tf], z_i, options);

    % Stato finale ottenuto
    Z_tf = Z(end,:)';
    
    % Condizione finale per le costate (lambda)
    lmb_tf = P*(Z_tf - z_f);

    % Integrazione delle equazioni aggiunte (costate) all’indietro nel tempo
    [Tlmb, lmb] = ode45(@(t,lmb) AdjointEq_A3(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,lmb_tf,m ,Cd, mu_r, epsilon), [tf t0], lmb_tf, options);

    % Ordinamento temporale crescente (per interpolazione)
    [Tlmb_sorted, idx] = sort(Tlmb);
    lmb_sorted = lmb(idx, :);

    % Interpolazione delle costate per farle corrispondere a Tz
    lmb1 = interp1(Tlmb_sorted, lmb_sorted(:,1), Tz);
    lmb2 = interp1(Tlmb_sorted, lmb_sorted(:,2), Tz);
    lmb3 = interp1(Tlmb_sorted, lmb_sorted(:,3), Tz);
    lmb4 = interp1(Tlmb_sorted, lmb_sorted(:,4), Tz);

    % Matrice delle costate
    LMB = [lmb1 lmb2 lmb3 lmb4]';
   
    % Estrazione degli stati
    Z1_Tz = Z(:,1);  Z2_Tz = Z(:,2);  Z3_Tz = Z(:,3);  Z4_Tz = Z(:,4);
    
    % Interpolazione del controllo sulla stessa griglia temporale
    u1_Tz = interp1(Tu, u(1,:), Tz, 'linear', 'extrap');
    u2_Tz = interp1(Tu, u(2,:), Tz, 'linear', 'extrap');

    % Funzione di costo per il vincolo morbido
    Soft_cost_fun = @(x,y) alpha.*exp((r^2 - (x-xc).^2 - (y-yc).^2)/sigma);

    % Lagrangiana istantanea (funzione di costo lungo la traiettoria)
    L_Tz = 0.5*(u1_Tz.*(w1*u1_Tz) + u2_Tz.*(w2*u2_Tz)) ...
           + Q*(Z4_Tz.^3) + Soft_cost_fun(Z1_Tz, Z2_Tz);

    % Funzione costo totale J (costo finale + integrale dei costi intermedi)
    J(ii,1) = 0.5*(Z(end,:).'-z_f)'*P*(Z(end,:).'-z_f) + trapz(Tz, L_Tz);

    % Derivata del gradiente del controllo dH/du (Hamiltoniana)
    [dH] = dHdu_new(u, Tu, lmb4, lmb3, Tz, m, R);

    % Norma del gradiente (criterio di convergenza)
    dH_norm = norm(dH, 'fro')

    % Controllo del criterio di arresto
    if dH_norm < eps
        break;  % convergenza raggiunta
    else 
        % Aggiornamento del controllo con passo "step"
        u_old = u;
        u = IterControl_new(dH, Tz, u_old, Tu, step);
        u(abs(u) < 1e-13) = 0;   % eliminazione numerica di valori molto piccoli
    end

end

%% ------------------------- Visualizzazione finale -----------------------
% Interpolazione della traiettoria finale sugli stessi punti del controllo
Z1 = interp1(Tz, Z(:,1), Tu);
Z2 = interp1(Tz, Z(:,2), Tu);

% Plot finale della traiettoria ottima
plot(x_cons, y_cons, 'k--', 'LineWidth', 1.5); 
hold on; grid on; grid minor;
scatter(pos_i(1), pos_i(2), 70, 'g', 'filled');
scatter(pos_f(1), pos_f(2), 70, 'r', 'filled'); 
plot(Z1, Z2, 'b', 'LineWidth', 2);
axis equal
legend('Vincolo morbido', 'Start', 'Goal', 'Traiettoria ottima');
title('Traiettoria finale dopo ottimizzazione');

%% plot J
% (plot opzionale per visualizzare la convergenza del costo)
% X = 1 : ii;
% semilogy(X, J(X,1))

%% LQR Regulator
% Re-sample the optimal trajectory 

% Time discretization in N = number of time intervals
N = 5000;

Ns = size(Z,2); % number of states
Nu = size(u,1); % number of controls

% Sampled time history
T = linspace(0,Tz(end),N);


% Interpolation/sampling of the state and control 
zk = zeros(N,Ns);
uk = zeros(N,Nu);
for ii = 1:Ns
    zk(:,ii) = interp1(Tz,Z(:,ii),T);
end
for ii = 1:Nu
    uk(:,ii) = interp1(Tu(1:end),u(ii,:),T);
end

% Option for ODE solver
VectTol = ones(Ns^2,1)*1e-5;
options = odeset('RelTol', 1e-5, 'AbsTol',VectTol);

% initial and final time
t0 = T(1);
tf = T(end);

% boundary conditions
p0 = P(1:end)';

% Integration of the matrix riccati equation
[Tp,PP] = ode23(@(t,p) DRE_A3(t,p,Q,R,zk,T,uk,T, Cd, m, mu_r, epsilon), flip(T), p0, options);

% From backward to forward dynamics (they are stored in the reversed order)
PP = flipud(PP);
Tp = flipud(Tp);

% Transformation Vector -> Matrix
PP_Matrix = zeros(Ns,Ns);

% Computation of the gain matrix in time, Uncontrolled stability matrix,
% Controlled stability matrix along the trajectory
K = zeros(N,Nu,Ns);
A = zeros(N,Ns,Ns);
Fc = zeros(N,Ns,Ns);

PolesUC = zeros(N,Ns);
PolesC = zeros(N,Ns);

Xp = zeros(N,Ns);

for ii = 1:N % Note: try also reshape.m
    % transformation vector -> matrix
    PP_Matrix(1:end) = PP(ii,:)';
    % control matrix G
    B = [0 0; 0 0; 0 1; 1/m 0];
    % Uncontrolled state stability matrix
    A(ii,:,:) = [0 0 -zk(ii,4)*sin(zk(ii,3))  cos(zk(ii,3));
         0 0  zk(ii,4)*cos(zk(ii,3))  sin(zk(ii,3));
         0 0  0           0;
         0 0  0       -2*Cd*zk(ii,4)/m- (mu_r*9.81)/(cosh(zk(ii,4)/epsilon)^2*epsilon)];
    % Uncontrolled system poles
    PolesUC(ii,:) = eig(squeeze(A(ii,:,:)));
    % Gain matrix C
    K(ii,:,:) = R^-1*B'*PP_Matrix; % R has to be non-null
    % Controlled state stability matrix
    Fc(ii,:,:) = squeeze(A(ii,:,:)) - B*squeeze(K(ii,:,:));
    % Controlled system poles
    PolesC(ii,:) = eig(squeeze(Fc(ii,:,:)));    
end


%% SIMULINK LQR
T = T';
IC = [0;0; 0; 0];
%K1 = squeeze(K);
K1 = K(:,1,:);
K1 = squeeze(K1);
% Wamp = 1; % amplitude of the disturbances

%% EKF

% observation matrix
C = [1,0,0,0;
    0,1,0,0];
Qk= eye(2);

Rk= [1,0;0,1];

Wamp = 2;
Wpamp = 2;
Wfreq = 10;      %[Hz]
Namp = 1;
out = sim('EKF_SImulink_Ass_3');
%% VANNO SISTEMATE LE LABEL
Xs = squeeze(out.Xs.data);
Xref = squeeze(out.Xref.data);
U = squeeze(out.U.data);
Uopt = squeeze(out.Uopt.data);
time = squeeze(out.time.data);

%% FIGURA 1: Stati (x, y, phi, v)
FigTag = figure;

% x position
subplot(2,2,1)
plot(time, Xref(1,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, Xs(:,1), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
ylabel('$x_{position}$ [m]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

% y position
subplot(2,2,2)
plot(time, Xref(2,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, Xs(:,2), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
ylabel('$y_{position}$ [m]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

% phi
subplot(2,2,3)
plot(time, Xref(3,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, Xs(:,3), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$\phi$ [rad]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

% v speed
subplot(2,2,4)
plot(time, Xref(4,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, Xs(:,4), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$v_{speed}$ [m/s]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

%% FIGURA 2: Ingressi di controllo (u1, u2)
FigTag = figure;

subplot(2,1,1)
plot(time, Uopt(1,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, U(1,:), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
ylabel('$u_{1}$ [N]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

subplot(2,1,2)
plot(time, Uopt(2,:), 'LineWidth', 1, 'LineStyle', '--', 'color', 'b'); hold on; grid on;
plot(time, U(2,:), 'LineWidth', 0.5, 'LineStyle', '-', 'color', 'r');
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$u_{2}$ [rad/s]', 'Interpreter', 'LaTex')
legend({'ref','sim'},'Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

%% FIGURA 3: Errori di stato (delta x, y, phi, v)
FigTag = figure;

subplot(2,2,1)
plot(time, Xs(:,1)' - Xref(1,:), 'LineWidth', 1, 'Color', 'b'); hold on; grid on;
ylabel('$\delta x$ [m]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

subplot(2,2,2)
plot(time, Xs(:,2)' - Xref(2,:), 'LineWidth', 1, 'Color', 'r'); hold on; grid on;
ylabel('$\delta y$ [m]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

subplot(2,2,3)
plot(time, Xs(:,3)' - Xref(3,:), 'LineWidth', 1, 'Color', 'g'); hold on; grid on;
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$\delta\phi$ [rad]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

subplot(2,2,4)
plot(time, Xs(:,4)' - Xref(4,:), 'LineWidth', 1, 'Color', 'm'); hold on; grid on;
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$\delta v$ [m/s]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

%% FIGURA 4: Errori di controllo (delta u1, u2)
FigTag = figure;

subplot(2,1,1)
plot(time, U(1,:) - Uopt(1,:), 'LineWidth', 1, 'Color', 'b'); hold on; grid on;
ylabel('$\delta u_{1}$ [N]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

subplot(2,1,2)
plot(time, U(2,:) - Uopt(2,:), 'LineWidth', 1, 'Color', 'r'); hold on; grid on;
xlabel('$t$ [s]', 'Interpreter', 'LaTex'); ylabel('$\delta u_{2}$ [rad/s]', 'Interpreter', 'LaTex')
legend('Error','Interpreter','Latex','Location','best')
ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'LaTex';

drawnow


%% --- PARAMETER ESTIMATION ---

% NON MALE 
% C = [1,0,0,0;
%     0,1,0,0];
% Qk= [1,0,0;0,1,0;0,0,1];
% Rk= [1,0;0,1];
% IC = [0;0; 0; 0];
% mu_r_guess = 1.2;
% Wamp = 0.5;
% Wfreq = 1000;      %[Hz]
% Namp = 0.5;


% observation matrix
C = [1,0,0,0;
    0,1,0,0];
Qk= [1,0,0;0,1,0;0,0,1];
Rk= [1,0;0,1];
IC = [0;0; 0; 0];
mu_r_guess = 1.2;
Wamp = 1;
Wfreq = 1000;      %[Hz]
Namp = 0.25;

out_param = sim("Param_Simulink_Ass_3.slx");

Xs = squeeze(out_param.Xs.data);
Xobs = squeeze(out_param.Xobs.data);
time = squeeze(out_param.time.data);

%% --- FIGURA A1: STIMA POSIZIONE (x, y) ---
% Confronto tra posizione reale (Xs) e stimata (Xobs)
FigTag = figure;

% Subplot X
ax = subplot(211);
h2 = plot(time, Xs(:,1), 'r'); hold on;
h2.LineWidth = 2;              % Linea più spessa per il vero stato
h2.Color = [1, 0, 0, 0.4];     % Rosso semi-trasparente
h1 = plot(time, Xobs(:,1), 'b--'); % Blu tratteggiato per la stima
h1.LineWidth = 1.2;
grid on
legend([h2,h1], {'$x_{true}$', '$\hat{x}_{est}$'}, 'Interpreter', 'LaTex', 'Location', 'best')
ylabel('$x$ [m]', 'Interpreter', 'LaTex')
ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% Subplot Y
ax = subplot(212);
h2 = plot(time, Xs(:,2), 'r'); hold on;
h2.LineWidth = 2;
h2.Color = [1, 0, 0, 0.4];
h1 = plot(time, Xobs(:,2), 'b--');
h1.LineWidth = 1.2;
grid on
legend([h2,h1], {'$y_{true}$', '$\hat{y}_{est}$'}, 'Interpreter', 'LaTex', 'Location', 'best')
xlabel('$t$ [s]', 'Interpreter', 'LaTex')
ylabel('$y$ [m]', 'Interpreter', 'LaTex')
ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% print(FigTag,'figA1_Pos.jpeg','-djpeg','-r600')

% --- FIGURA A2: STIMA DINAMICA (phi, v) ---
% Confronto tra angoli e velocità
FigTag = figure;

% Subplot Phi
ax = subplot(211);
h2 = plot(time, Xs(:,3), 'r'); hold on;
h2.LineWidth = 2;
h2.Color = [1, 0, 0, 0.4];
h1 = plot(time, Xobs(:,3), 'b--');
h1.LineWidth = 1.2;
grid on
legend([h2,h1], {'$\phi_{true}$', '$\hat{\phi}_{est}$'}, 'Interpreter', 'LaTex', 'Location', 'best')
ylabel('$\phi$ [rad]', 'Interpreter', 'LaTex')
ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% Subplot Velocità
ax = subplot(212);
h2 = plot(time, Xs(:,4), 'r'); hold on;
h2.LineWidth = 2;
h2.Color = [1, 0, 0, 0.4];
h1 = plot(time, Xobs(:,4), 'b--');
h1.LineWidth = 1.2;
grid on
legend([h2,h1], {'$v_{true}$', '$\hat{v}_{est}$'}, 'Interpreter', 'LaTex', 'Location', 'best')
xlabel('$t$ [s]', 'Interpreter', 'LaTex')
ylabel('$v$ [m/s]', 'Interpreter', 'LaTex')
ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% print(FigTag,'figA2_Dyn.jpeg','-djpeg','-r600')

% --- FIGURA A3: STIMA PARAMETRO (Attrito) ---
% Questo è il plot più importante per il tuo obiettivo
FigTag = figure;
ax = axes;

% Plot del vero parametro (che potrebbe variare o essere costante)
h2 = yline(mu_r, 'r'); hold on;
h2.LineWidth = 2.5;
h2.Color = [1, 0, 0, 0.4]; % Rosso trasparente per il valore vero

% Plot della stima del parametro
h1 = plot(time, Xobs(:,5), 'b'); 
h1.LineWidth = 1.5;

grid on
xlabel('$t$ [s]', 'Interpreter', 'LaTex')
ylabel('Friction $\mu$', 'Interpreter', 'LaTex')
title('Parameter Estimation: Friction Coefficient', 'Interpreter', 'LaTex')
legend([h2,h1], {'$\mu_{true}$', '$\hat{\mu}_{est}$'}, 'Interpreter', 'LaTex', 'Location', 'best')

% Zoom automatico per vedere meglio la convergenza
ylim([min(Xobs(:,5))*0.8, max(Xobs(:,5))*1.2])

ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% print(FigTag,'figA3_Param.jpeg','-djpeg','-r600')

% --- FIGURA A4: ERRORE DI STIMA PARAMETRO ---
FigTag = figure;
ax = axes;

% Calcolo errore
param_error = (mu_r*ones(length(Xobs),1) - Xobs(:,5))/mu_r *100;

h1 = plot(time, param_error, 'k'); % Nero per l'errore
h1.LineWidth = 1.5;
grid on; hold on;

% Linea dello zero per riferimento
yline(0, 'r--', 'LineWidth', 1);

xlabel('$t$ [s]', 'Interpreter', 'LaTex')
ylabel('$ \%  Relative  \mu$ error', 'Interpreter', 'LaTex')
legend({'$(\mu_{true} - \hat{\mu}_{est})$'}, 'Interpreter', 'LaTex', 'Location', 'best')
title('Parameter Estimation Convergence Error', 'Interpreter', 'LaTex')

ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% print(FigTag,'figA4_Err.jpeg','-djpeg','-r600')