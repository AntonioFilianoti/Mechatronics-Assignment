function dH = dHdu(u, Tu, LMB, t, m, R)
% dHdu - Calcola il gradiente della Hamiltoniana rispetto al controllo
%
% INPUT:
%   u     - matrice controlli [2xN]
%   Tu    - griglia temporale
%   LMB   - matrice costate [4xN]
%   m,R   - parametri dinamici e peso controllo
%
% OUTPUT:
%   dH    - gradiente Hamiltoniana rispetto a u [2xN]

 u1 = interp1(Tu, u(1,:), t);
 u2 = interp1(Tu, u(2,:), t);

    L_u = R*u;  % derivata della parte quadratica del costo
    B  = [0 0;
          0 0;
          0 1;
          1/m 0];

    dH = L_u + B' * LMB;
end

