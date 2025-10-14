function dz = stateEq(t,z,u,Tu,m,Cd)
% stateEq - Equazioni di stato del sistema dinamico
%
% INPUT:
%   t   - tempo attuale [s]
%   z   - vettore di stato [x; y; theta; v]
%   u   - matrice di controlli [2 x N]
%   Tu  - vettore dei tempi associati a u
%   m   - massa
%   Cd  - coefficiente di resistenza
%
% OUTPUT:
%   dz  - derivate di stato [4x1]

    u1 = interp1(Tu, u(1,:), t);
    u2 = interp1(Tu, u(2,:), t);

    dz = zeros(4,1);
    dz(1) = z(4)*cos(z(3));          % x_dot
    dz(2) = z(4)*sin(z(3));          % y_dot
    dz(3) = u2;                      % theta_dot
    dz(4) = (1/m)*(u1 - Cd*z(4)^2);  % v_dot
end
