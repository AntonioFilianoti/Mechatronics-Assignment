function [dz] = stateEq_A3(t,z,u,Tu, m, Cd, mu_r, epsilon)

%% ========================================================================
% Funzione: stateEq_new
% ------------------------------------------------------------------------
% Descrizione:
%   Definisce le equazioni di stato del sistema dinamico non lineare.
%   Lo stato è composto da:
%       z = [x; y; θ; v]
%   dove:
%       x, y  → posizione
%       θ     → orientamento
%       v     → velocità
%   I controlli sono:
%       u1 → forza di trazione (accelerazione lungo la direzione θ)
%       u2 → controllo di rotazione (variazione dell'angolo θ)
%   Il modello include una resistenza aerodinamica proporzionale a v².
%
%   Equazioni:
%       ẋ = v*cos(θ)
%       ẏ = v*sin(θ)
%       θ̇ = u2
%       v̇ = (1/m)*(u1 - Cd*v²)
%
% Input:
%   t   → tempo corrente
%   z   → stato corrente [x; y; θ; v]
%   u   → controllo definito sulla griglia temporale Tu
%   Tu  → vettore temporale del controllo
%   m   → massa
%   Cd  → coefficiente di resistenza aerodinamica
%
% Output:
%   dz  → derivata dello stato (vettore 4x1)
%% ========================================================================

    dz = zeros(4,1);                     % inizializza derivata stato

    % Interpolazione dei controlli u1(t) e u2(t) al tempo corrente t
    u1 = interp1(Tu,u(1,:),t);
    u2 = interp1(Tu,u(2,:),t);

    % Equazioni dinamiche del sistema
    dz(1) = z(4)*cos(z(3));              % ẋ = v * cos(θ)
    dz(2) = z(4)*sin(z(3));              % ẏ = v * sin(θ)
    dz(3) = u2;                          % θ̇ = u2 (velocità angolare)
    dz(4) = 1/m*(u1 - Cd*z(4)^2 - mu_r*m*9.81*tanh(z(4)/epsilon));        % v̇ = (u1 - Cd*v²)/m
end

