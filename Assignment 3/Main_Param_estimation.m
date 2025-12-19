clc
clear
close all


%% --- PARAMETER ESTIMATION ---
mu_r = 1;
mu_r_guess = 3;

% observation matrix
C = [1,0,0,0;
    0,1,0,0];
Wamp = 0;
WAamp = 2;
Wfreq = 1000;      %[Hz]
Namp = 0.01;
n = 0.01;
Qk= diag([Wamp Wamp WAamp]);
Rk= [n,0;0,n];
IC = [0;0; 0; 0];
tf = 100;


out_param = sim("Param_Sinewave_Simulink_Ass_3.slx");

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
ylabel('$ \% - Relative - \mu$ error', 'Interpreter', 'LaTex')
legend({'$(\mu_{true} - \hat{\mu}_{est})$'}, 'Interpreter', 'LaTex', 'Location', 'best')
title('Parameter Estimation Convergence Error', 'Interpreter', 'LaTex')

ax.FontSize = 16;
ax.TickLabelInterpreter = 'LaTex';

% print(FigTag,'figA4_Err.jpeg','-djpeg','-r600')