function [dH] = dHdu_new(u, Tu, lmb4, lmb3, Tz, m, R)

%% ========================================================================
% Funzione: dHdu_new
% ------------------------------------------------------------------------
% Descrizione:
%   Calcola la derivata dell’Hamiltoniana rispetto ai controlli u.
%
%   L’Hamiltoniana H è data da:
%       H = L + λᵀ f(z, u)
%   dove:
%       L → costo istantaneo
%       λ → vettore delle costate (aggiunte)
%       f → dinamica del sistema
%
%   La derivata dH/du rappresenta il gradiente del costo rispetto ai controlli.
%
% Input:
%   u       → controllo sulla griglia Tu
%   Tu      → tempi del controllo
%   lmb4    → quarta costata λ₄ (associata alla velocità)
%   lmb3    → terza costata λ₃ (associata all’orientamento)
%   Tz      → tempi della simulazione degli stati
%   m       → massa
%   R       → matrice dei pesi sui controlli
%
% Output:
%   dH      → derivata dell’Hamiltoniana rispetto a u = [u1; u2]
%% ========================================================================

    % Interpolazione dei controlli u1(t) e u2(t) alla griglia Tz
    u1 = interp1(Tu,u(1,:),Tz);
    u2 = interp1(Tu,u(2,:),Tz);
    u = [u1 u2]';                        % ricostruisce matrice dei controlli

    % Matrice B del sistema (relazione tra controllo e derivata dello stato)
    B = [0 0; 0 0; 0 1; 1/m 0];

    % Calcolo del gradiente dell’Hamiltoniana:
    % dH/du₁ = ∂L/∂u₁ + (∂f/∂u₁)ᵀ λ = R₁₁*u₁ + λ₄/m
    % dH/du₂ = ∂L/∂u₂ + (∂f/∂u₂)ᵀ λ = R₂₂*u₂ + λ₃
    dH = [u1*R(1,1) + (lmb4/m), lmb3 + u2*R(2,2)]';
end