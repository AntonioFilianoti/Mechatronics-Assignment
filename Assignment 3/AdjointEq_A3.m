function dlmb = AdjointEq_new(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m, lmb_tf )
%% ========================================================================
% Funzione: AdjointEq_new
% ------------------------------------------------------------------------
% Descrizione:
%   Calcola le **equazioni aggiunte** (o equazioni delle costate) λ̇(t)
%   per un sistema dinamico non lineare nel contesto del controllo ottimo.
%
%   Le equazioni aggiunte derivano dalle condizioni di ottimalità di Pontryagin:
%
%       λ̇ = - ( ∂H/∂z )ᵗ
%
%   dove H è l’Hamiltoniana:
%       H = L(z,u) + λᵀ f(z,u)
%
%   Qui:
%       - L(z,u) è la funzione di costo istantaneo
%       - f(z,u) rappresenta la dinamica del sistema
%       - λ è il vettore delle costate [λ₁, λ₂, λ₃, λ₄]
%
% Input:
%   t        → tempo corrente
%   lmb      → vettore delle costate λ(t)
%   Z        → matrice degli stati simulati nel tempo [x, y, θ, v]
%   Tz       → vettore dei tempi corrispondenti agli stati Z
%   Q        → peso sul termine di penalizzazione della velocità
%   sigma    → parametro di “morbidezza” del vincolo morbido
%   alpha    → peso del vincolo morbido
%   xc, yc   → coordinate del centro del vincolo morbido
%   r        → raggio dell’area vincolata
%   Cd       → coefficiente di resistenza aerodinamica
%   m        → massa
%   lmb_tf   → valore finale (condizione al contorno) di λ(t_f)
%
% Output:
%   dlmb     → derivata temporale delle costate (vettore 4x1)
%% ========================================================================

    dlmb = zeros(4,1);  % inizializza il vettore derivata costate (λ̇)

    % --------------------------------------------------------------------
    % Interpolazione degli stati z(t) al tempo corrente t
    % --------------------------------------------------------------------
    z1 = interp1(Tz,Z(:,1), t);   % posizione x
    z2 = interp1(Tz,Z(:,2), t);   % posizione y
    z3 = interp1(Tz,Z(:,3), t);   % orientamento θ
    z4 = interp1(Tz,Z(:,4), t);   % velocità v

    % (Eventuali interpolazioni dei controlli u1, u2 commentate)
    % u1 = interp1(Tu,u(1,:),t);
    % u2 = interp1(Tu,u(2,:),t);

    % --------------------------------------------------------------------
    % Matrice Jacobiana A = ∂f/∂z
    %   Rappresenta la derivata della dinamica rispetto agli stati.
    %   Serve per il termine Aᵗ·λ nelle equazioni aggiunte.
    % --------------------------------------------------------------------
    A = [0 0 -z4*sin(z3)  cos(z3);
         0 0  z4*cos(z3)  sin(z3);
         0 0  0           0;
         0 0  0       -2*Cd*z4/m];

    % --------------------------------------------------------------------
    % Calcolo del gradiente del vincolo morbido rispetto a (x,y)
    %   Il vincolo è modellato con una funzione esponenziale che cresce
    %   vicino alla zona proibita (soft constraint).
    %
    %   Vincolo: alpha * exp((r² - (x - xc)² - (y - yc)²) / sigma)
    % --------------------------------------------------------------------
    expo = (r^2 - (z1-xc)^2 - (z2-yc)^2) / sigma;
    exp_term = alpha * exp(expo);

    % Derivate parziali rispetto a x e y
    a = -(2*(z1 - xc)/sigma) * exp_term;  % ∂L/∂x
    b = -(2*(z2 - yc)/sigma) * exp_term;  % ∂L/∂y

    % --------------------------------------------------------------------
    % Gradiente del Lagrangiano rispetto allo stato z
    %   Lz = [∂L/∂x, ∂L/∂y, ∂L/∂θ, ∂L/∂v]
    % --------------------------------------------------------------------
    Lz = [a, b, 0, 3*Q*z4^2];  % il termine su v deriva da Q*v³ nel costo

    % --------------------------------------------------------------------
    % Equazioni aggiunte (costate):
    %   λ̇ = - [ ∂L/∂z + (∂f/∂z)ᵗ λ ]
    % --------------------------------------------------------------------
    dlmb = -(Lz + (lmb' * A))';

    % (alternativa equivalente, lasciata come commento nel codice originale)
    % dlmb = -Lz' - A * lmb;
end