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
m   = 1;          % massa [kg]
Cd  = 0.15;       % coefficiente di resistenza
t0  = 0;          % tempo iniziale [s]
tf  = 1.6;         % tempo finale [s]

% Stati: [x; y; theta; v]
pos_i = [0; 0; 0];  
pos_f = [1;1; pi/3];

vel_i = 0;                 % velocità iniziale
vel_f = 0;                 % velocità finale

z_i = [pos_i; vel_i];
z_f = [pos_f; vel_f];

% max state vector used for normalization/scaling
xmax_vec = [pos_f(1:2); 2*pi/3;1.5];
% max control vector used for normalization/scaling
umax_vec = [1, 1]';
%% ----------------------- Pesi e vincoli del costo -----------------------
w1 = 0.2;
w2 = 0.2;
%final state weight
P = diag(1./xmax_vec.^2);
P = P.*diag([180 180 20 0]);
%optimal control weight
R = diag(1./umax_vec.^2);
R = R.*diag([w1 w2]);

Q = 0.05;                       % peso su penalizzazione velocità

alpha = 2;                    % peso vincolo morbido (soft constraint)
sigma = 0.01;                     % parametro di “softness”
xc = 0.65;
yc = 0.65;
r = 0.2;
            

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
eps  = 5e-3;                           % tolleranza di arresto
u = [ones(1, Nsegment); zeros(1, Nsegment)]; % controllo iniziale

%% ------------------------- Procedura Iterativa --------------------------
for ii = 1:Nmax

    [Tz, Z] = ode45(@(t,z) stateEq_new(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);

    % Z(abs(Z) < 1e-15) = 0;

    Z_tf = Z(end,:)';
    
   lmb_tf = P*(Z_tf-z_f);

   [Tlmb, lmb] = ode45(@(t,lmb) AdjointEq_new(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m,lmb_tf  ), [tf t0], lmb_tf, options);
  [Tlmb_sorted, idx] = sort(Tlmb);   % increasing times
    lmb_sorted = lmb(idx, :);

    % Now interpolate using increasing time vector
    lmb1 = interp1(Tlmb_sorted, lmb_sorted(:,1), Tz);
    lmb2 = interp1(Tlmb_sorted, lmb_sorted(:,2), Tz);
    lmb3 = interp1(Tlmb_sorted, lmb_sorted(:,3), Tz);
    lmb4 = interp1(Tlmb_sorted, lmb_sorted(:,4), Tz);


   LMB = [lmb1 lmb2 lmb3 lmb4]';
   

%-----cost function----------
Z1_Tz = Z(:,1);  Z2_Tz = Z(:,2);  Z3_Tz = Z(:,3);  Z4_Tz = Z(:,4);
% u(t) sulla griglia Tz
u1_Tz = interp1(Tu, u(1,:), Tz, 'linear', 'extrap');
u2_Tz = interp1(Tu, u(2,:), Tz, 'linear', 'extrap');

Soft_cost_fun = @(x,y) alpha.*exp((r^2 - (x-xc).^2 - (y-yc).^2)/sigma);

L_Tz = 0.5*(u1_Tz.*(w1*u1_Tz) + u2_Tz.*(w2*u2_Tz)) ...
       + Q*(Z4_Tz.^3) + Soft_cost_fun(Z1_Tz, Z2_Tz);

J(ii,1) = 0.5*(Z(end,:).'-z_f)'*P*(Z(end,:).'-z_f) + trapz(Tz, L_Tz);

% ---dH----
  [dH ]= dHdu_new(u, Tu,lmb4, lmb3,Tz, m, R);
  %dH_1_norm = norm(dH_1, 'fro')
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

%% plot J

%X = 1 : ii;
%semilogy(X, J(X,1))
