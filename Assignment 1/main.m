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
global m Cd
m   = 1;          % massa del corpo [kg]
Cd  = 0.15;       % coefficiente di resistenza aerodinamica (attrito quadratico)
t0  = 0;          % tempo iniziale [s]
tf  = 1.6;        % tempo finale [s]

% Stati del sistema: [x; y; theta; v]
pos_i = [0; 0; 0];         % posizione e orientamento iniziali (x, y, theta)
pos_f = [1; 1; pi/3];      % posizione e orientamento finali desiderati

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
xlim([-2 12]); ylim([-2 12]);
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
eps  = 5e-3;                            % soglia di convergenza
u = [ones(1, Nsegment); zeros(1, Nsegment)];  % controllo iniziale (u1=1, u2=0)

%% ------------------------- Procedura Iterativa --------------------------
for ii = 1:Nmax

    % Integrazione delle equazioni di stato in avanti nel tempo
    [Tz, Z] = ode45(@(t,z) stateEq_new(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);

    % Stato finale ottenuto
    Z_tf = Z(end,:)';
    
    % Condizione finale per le costate (lambda)
    lmb_tf = P*(Z_tf - z_f);

    % Integrazione delle equazioni aggiunte (costate) all’indietro nel tempo
    [Tlmb, lmb] = ode45(@(t,lmb) AdjointEq_new(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m,lmb_tf), [tf t0], lmb_tf, options);

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
[Tp,PP] = ode23(@(t,p) DRE(t,p,Q,R,zk,T,uk,T), flip(T), p0, options);

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
         0 0  0       -2*Cd*zk(ii,4)/m];
    % Uncontrolled system poles
    PolesUC(ii,:) = eig(squeeze(A(ii,:,:)));
    % Gain matrix C
    K(ii,:,:) = R^-1*B'*PP_Matrix; % R has to be non-null
    % Controlled state stability matrix
    Fc(ii,:,:) = squeeze(A(ii,:,:)) - B*squeeze(K(ii,:,:));
    % Controlled system poles
    PolesC(ii,:) = eig(squeeze(Fc(ii,:,:)));    
end

FigTag = figure;
FigPos = FigTag.Position;
FigTag.Position = [FigPos(1)-FigPos(3)*0.25,FigPos(2)-FigPos(4)*0.2,FigPos(3)*1.5,FigPos(4)*1.2];
filename = 'Trajectory.gif';
PlotSpeed = 20;
for ii = 1:N/PlotSpeed
    
    clf(FigTag)
    
    % Temporal evolution - uncontrolled system
    ax = subplot(2,2,1); hold on; grid on;
    plot(real(PolesUC(ii*PlotSpeed,:)),imag(PolesUC(ii*PlotSpeed,:)),'x',...
                'LineWidth',2,'Color','k'); hold on; box on;
    eigRe = real(PolesUC(1:ii*PlotSpeed,:));
    eigIm = imag(PolesUC(1:ii*PlotSpeed,:));
    Tmatrix = Tp(1:ii*PlotSpeed)*ones(1,size(eigRe,2));
    scatter(eigRe(:),eigIm(:),2,Tmatrix(:),'fill'); hold on; box on;
    colormap jet; caxis([0,Tp(end)]);
    title('Uncontrolled system','Interpreter','LaTex')
    xlabel('$\Re\{\lambda\}$ [Rad/s]','Interpreter','LaTex');
    ylabel('$\Im\{\lambda\}$ [Rad/s]','Interpreter','LaTex');
    ax.FontSize = 14;
    ax.TickLabelInterpreter = 'LaTex';
    ylim([-25,25])
    xlim([-50,10])
    
    % Temporal evolution - controlled system
    ax = subplot(2,2,2); hold on; grid on;
    plot(real(PolesC(ii*PlotSpeed,:)),imag(PolesC(ii*PlotSpeed,:)),'x',...
                'LineWidth',2,'Color','k'); hold on; box on;
    eigRe = real(PolesC(1:ii*PlotSpeed,:));
    eigIm = imag(PolesC(1:ii*PlotSpeed,:));
    Tmatrix = Tp(1:ii*PlotSpeed)*ones(1,size(eigRe,2));
    scatter(eigRe(:),eigIm(:),2,Tmatrix(:),'fill'); hold on; box on;
    colormap jet;  caxis([0,Tp(end)]);
    title('Controlled system','Interpreter','LaTex')
    xlabel('$\Re\{\lambda\}$ [Rad/s]','Interpreter','LaTex');
    ylabel('$\Im\{\lambda\}$ [Rad/s]','Interpreter','LaTex');
    ax.FontSize = 14;
    ax.TickLabelInterpreter = 'LaTex';
    ylim([-25,25])
    xlim([-50,10])

    
    % plot of the trajectory
    ax = subplot(2,2,[3,4]);
    hold on; grid on;
    title('Trajectory','Interpreter','LaTex')
    scatter(T(1:ii*PlotSpeed),zk(1:ii*PlotSpeed,1),2,T(1:ii*PlotSpeed),'fill');
    scatter(T(1:ii*PlotSpeed),zk(1:ii*PlotSpeed,2),2,T(1:ii*PlotSpeed),'fill');
    scatter(T(1:ii*PlotSpeed),zk(1:ii*PlotSpeed,3),2,T(1:ii*PlotSpeed),'fill');
    scatter(T(1:ii*PlotSpeed),zk(1:ii*PlotSpeed,4),2,T(1:ii*PlotSpeed),'fill');
    colormap jet; caxis([0,Tp(end)]);
    xlabel('$t$ [s]','Interpreter','LaTex');
    ylabel(',\theta$','Interpreter','LaTex');
    plot(T(ii*PlotSpeed),zk(ii*PlotSpeed,1),'o','LineWidth',1,'color','k');
    plot(T(ii*PlotSpeed),zk(ii*PlotSpeed,2),'o','LineWidth',1,'color','k');
    ax.FontSize = 14;
    ax.TickLabelInterpreter = 'LaTex';
    hold on; box on;
    xlim([0,T(end)])
    ylim([-max(max(abs(zk))),max(max(abs(zk)))])
    
    set(gcf,'Color','w')
%     axis equal
    drawnow    
    
    % % Save Animation
    % frame = getframe(FigTag);
    % im = frame2im(frame); 
    % [imind,cm] = rgb2ind(im,256); 
    % if ii==1
    %  imwrite(imind,cm,filename,'gif','DelayTime',0.001, 'Loopcount',inf); 
    % else 
    %  imwrite(imind,cm,filename,'gif','DelayTime',0.001,'WriteMode','append'); 
    % end         
end
    