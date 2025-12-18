function u = IterControl_new(dH, Tz, u_old, Tu, step)
%% ========================================================================
% Funzione: IterControl_new
% ------------------------------------------------------------------------
% Descrizione:
%   Aggiorna il controllo in base al gradiente del costo (dH/du).
%   È il passo iterativo principale del metodo del gradiente.
%
%   Formula:
%       u_new = u_old - step * dH
%
%   Dove "step" è il passo di discesa.
%
% Input:
%   dH      → derivata parziale dell’Hamiltoniana rispetto ai controlli
%   Tz      → tempi corrispondenti al vettore dH
%   u_old   → controllo precedente (matrice 2xN)
%   Tu      → griglia temporale del controllo
%   step    → passo di aggiornamento
%
% Output:
%   u       → nuovo controllo aggiornato
%% ========================================================================

    % Interpola il gradiente dH alla griglia temporale del controllo
    dH1 = interp1(Tz,dH(1,:), Tu);
    dH2 = interp1(Tz,dH(2,:), Tu);

    dH = [dH1; dH2];                     % ricompone la matrice del gradiente

    % Aggiornamento del controllo secondo la discesa del gradiente
    u = u_old - step*dH;
end