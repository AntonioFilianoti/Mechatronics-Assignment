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
m   = 1;          % massa [kg]
Cd  = 0.35;       % coefficiente di resistenza
t0  = 0;          % tempo iniziale [s]
tf  = 10;         % tempo finale [s]

% Stati: [x; y; theta; v]
pos_i = [0; 0; pi/4];      % posizione e orientamento iniziali
pos_f = [10; 10; 0];       % posizione e orientamento finali
vel_i = 0;                 % velocità iniziale
vel_f = 1;                 % velocità finale

z_i = [pos_i; vel_i];
z_f = [pos_f; vel_f];

%% ----------------------- Pesi e vincoli del costo -----------------------
R = diag([0.0001, 0.0001]);             % peso sul controllo
Q = 0.05;                       % peso su penalizzazione velocità
P = diag([1000, 1000, 1, 1]);     % peso su stato finale
alpha = 10;                    % peso vincolo morbido (soft constraint)
sigma = 10;                     % parametro di “softness”
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
Nsegment = 500;                       % numero di passi temporali
Tu = linspace(t0, tf, Nsegment);       % vettore temporale
options = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% ----------------------- Parametri dell’iterazione ----------------------
Nmax = 5000;                           % massimo numero iterazioni
step = 5e-5;                           % passo di aggiornamento controllo
eps  = 1e-2;                           % tolleranza di arresto
u = [ones(1, Nsegment); zeros(1, Nsegment)]; % controllo iniziale

%% ------------------------- Procedura Iterativa --------------------------
for ii = 1:Nmax

    % 1) Integrazione avanti degli stati
    tic
    [Tz, Z] = ode45(@(t,z) stateEq(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);
    toc

    % 2) Integrazione indietro delle aggiunte
    Z_tf = Z(end,:)';
    lmb_tf = P*(Z_tf - z_f);
    [Tlmb, lmb] = ode45(@(t,l) AdjointEq(t,l,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m), [t0 tf], lmb_tf, options);

    % Ribalta i risultati per allineare i tempi
    lmb = sort(lmb);
    Tlmb = sort(Tlmb);

    % Interpola le aggiunte per allinearle ai tempi di stato
    LMB = interp1(Tlmb, lmb, Tu, 'linear', 'extrap')';
    
    % 3) Calcolo del gradiente della Hamiltoniana
    dH = dHdu(u, Tu, LMB, Tz, m, R);
    H_Norm = norm(dH, 'fro');

    % Evita divergenza numerica
    if H_Norm >= 1e10
        H_Norm = 1e9;
        dH = dH ./ 1e6;
    end

    % 4) Aggiornamento controllo (gradiente discendente)
    u = u - step * dH;

    % Condizione di uscita opzionale
    if H_Norm < eps
        disp(['Convergenza raggiunta in ', num2str(ii), ' iterazioni.']);
        break;
    end

    a=0;
end

%% ------------------------- Visualizzazione finale -----------------------
Z1 = interp1(Tz, Z(:,1), Tu);
Z2 = interp1(Tz, Z(:,2), Tu);

plot(Z1, Z2, 'b', 'LineWidth', 2);
legend('Vincolo morbido', 'Start', 'Goal', 'Traiettoria ottima');
title('Traiettoria finale dopo ottimizzazione');



