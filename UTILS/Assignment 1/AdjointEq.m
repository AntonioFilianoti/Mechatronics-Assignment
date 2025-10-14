function dlmb = AdjointEq(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m)
% AdjointEq - Calcola le equazioni aggiunte (costate)
%
% INPUT:
%   t    - tempo attuale
%   lmb  - vettore costate [4x1]
%   Z,Tz - traiettoria di stato e tempi
%   Q, sigma, alpha, xc, yc, r - parametri del costo e vincolo
%   Cd, m - parametri dinamici
%
% OUTPUT:
%   dlmb - derivate delle costate [4x1]

    % Interpolazione degli stati
    z = interp1(Tz, Z, t);
    x=z(1); y=z(2); th=z(3); v=z(4);

    % Matrice Jacobiana della dinamica
    A = [0 0 -v*sin(th) cos(th);
         0 0  v*cos(th) sin(th);
         0 0  0          0;
         0 0  0       -2*Cd*v/m];

    % Gradiente del vincolo morbido
    expo = (r^2 - (x-xc)^2 - (y-yc)^2)/sigma;
    expo = min(expo, 1000); % previene overflow
    exp_term = exp(expo);
    a = -alpha*(2*(x-xc)/sigma)*exp_term;
    b = -alpha*(2*(y-yc)/sigma)*exp_term;

    % Gradiente del Lagrangiano rispetto a x
    Lx = [a b 0 3*Q*v^2];
    dlmb = -(Lx + (lmb'*A))';
end
