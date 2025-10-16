clear
clc


%% ========================================================================
%  Optimal Control with Adjoint Equations (Nonlinear System)
%  ------------------------------------------------------------------------
%  Questo script risolve un problema di controllo ottimo per un sistema 
%  dinamico non lineare (con attrito quadratico e vincolo morbido di traiettoria)
%  utilizzando un approccio iterativo basato su equazioni aggiunte (costate).
%
%  ------------------------------------------------------------------------
%% ========================================================================

clear; close all; clc;

%% ------------------------ Parametri del modello -------------------------
m   = 100;          % massa [kg]
Cd  = 0.35;       % coefficiente di resistenza
t0  = 0;          % tempo iniziale [s]
tf  = 25*4;         % tempo finale [s]

% Stati: [x; y; theta; v]
pos_i = [0; 0; pi/4];  
pos_f = [10; 10; pi/4];
% pos_i = [0; 0;0];      % posizione e orientamento iniziali
% pos_f = [10; 0; 0];       % posizione e orientamento finali
vel_i = 0;                 % velocità iniziale
vel_f = 0;                 % velocità finale

z_i = [pos_i; vel_i];
z_f = [pos_f; vel_f];

%% ----------------------- Pesi e vincoli del costo -----------------------
R = diag([1, 1]);             % peso sul controllo
Q = 0;                       % peso su penalizzazione velocità
P = diag([200, 200, 1, 100]);     % peso su stato finale
alpha = 0;                    % peso vincolo morbido (soft constraint)
sigma = 0.1;                     % parametro di “softness”
r = 1;                          % raggio vincolo circolare
xc = 5; yc = 5;                 % centro vincolo

%% -------------------------- Plot del vincolo ----------------------------
theta = linspace(0, 2*pi, 400);
x_cons = xc + r .* cos(theta);
y_cons = yc + r .* sin(theta);

figure('Position', [100 100 900 700]);
plot(x_cons, y_cons, 'k--', 'LineWidth', 1.5); hold on;
scatter(pos_i(1), pos_i(2), 70, 'g', 'filled');
scatter(pos_f(1), pos_f(2), 70, 'r', 'filled');
xlabel('x'); ylabel('y');
axis equal;
grid on; grid minor;
xlim([-2 12]); ylim([-2 12]);
legend('Vincolo morbido', 'Start', 'Goal');
title('Problema di controllo ottimo con vincolo morbido');

%% --------------------------- Discretizzazione ---------------------------
Nsegment = 2000;                       % numero di passi temporali
Tu = linspace(t0, tf, Nsegment);       % vettore temporale
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
options_state = odeset('RelTol',1e-12,'AbsTol',1e-14);

%% ----------------------- Parametri dell’iterazione ----------------------
Nmax = 5*10000;                           % massimo numero iterazioni
step = 5e-3;                           % passo di aggiornamento controllo
eps  = 1e-2;                           % tolleranza di arresto
u = [ones(1, Nsegment); zeros(1, Nsegment)]; % controllo iniziale

%% ------------------------- Procedura Iterativa --------------------------
for ii = 1%:Nmax

    [Tz, Z] = ode45(@(t,z) stateEq_new(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);

    % Z(abs(Z) < 1e-15) = 0;

    Z_tf = Z(end,:).';
    
   lmb_tf = P*(Z_tf-z_f);

   [Tlmb, lmb] = ode45(@(t,lmb) AdjointEq_new(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m,lmb_tf  ), [tf t0], lmb_tf, options);
   Tlmb = flipud(Tlmb); 
   lmb = flipud(lmb);

   lmb1 = interp1(Tlmb, lmb(:,1), Tz);
   lmb2 = interp1(Tlmb, lmb(:,2), Tz);
   lmb3 = interp1(Tlmb, lmb(:,3), Tz);
   lmb4 = interp1(Tlmb, lmb(:,4), Tz);

   LMB = [lmb1 lmb2 lmb3 lmb4]';
   %LMB(abs(LMB) < 1e-11) = 0;
   dH = dHdu_new(u, Tu, LMB, Tz, m, R);
   dH_norm = norm(dH, 'fro')

   if dH_norm < eps
       break;
   else 
       u_old = u;
       u = IterControl_new(dH, Tz, u_old, Tu, step);
       u(abs(u) < 1e-13) = 0;
   end



end

%% ------------------------- Visualizzazione finale -----------------------
Z1 = interp1(Tz, Z(:,1), Tu);
Z2 = interp1(Tz, Z(:,2), Tu);


plot(x_cons, y_cons, 'k--', 'LineWidth', 1.5); 
hold on; grid on; grid minor;
scatter(pos_i(1), pos_i(2), 70, 'g', 'filled');
scatter(pos_f(1), pos_f(2), 70, 'r', 'filled'); 
plot(Z1, Z2, 'b', 'LineWidth', 2);
legend('Vincolo morbido', 'Start', 'Goal', 'Traiettoria ottima');
title('Traiettoria finale dopo ottimizzazione');