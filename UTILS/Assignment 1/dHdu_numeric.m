function dH = dHdu_numeric(u, Tu, LMB, Tz, m, R)
% dHdu_numeric: stima numerica di dH/du (Hamiltonian gradient wrt controls)
%
% Inputs:
%   u    -> 2xN control vector
%   Tu   -> time grid (1xN)
%   LMB  -> 4xN matrix of adjoint vectors (interpolated)
%   Tz   -> time grid of states (usato solo per consistenza)
%   m, R -> parametri dinamica / costi
%
% Output:
%   dH   -> 2xN matrix of numerical partials dH/du

    eps_fd = 1e-6;          % passo differenze finite
    dH = zeros(size(u));    % prealloca

    for k = 1:length(Tu)
        uk = u(:,k);
        lam = LMB(k,:);
        % Hamiltonian locale
        H_fun = @(uc) Hamiltonian_local(uc, lam, m, R);

        % differenze finite centrali
        for i = 1:2
            du = zeros(2,1);
            du(i) = eps_fd;
            Hp = H_fun(uk + du);
            Hm = H_fun(uk - du);
            dH(i,k) = (Hp - Hm)/(2*eps_fd);
        end
    end
end


function H = Hamiltonian_local(u, lmb, m, R)
% Hamiltoniana locale (solo parte di controllo)
% Assumiamo che la dinamica dipenda dai comandi come:
%   dot{z} = [v*cos(theta);
%              v*sin(theta);
%              u2/m;
%              u1/m]
% (aggiorna secondo la tua "stateEq")

    % controllo lineare in cost
    H = 0.5*(u.'*R*u) + lmb(4)*(u(1)/m) + lmb(3)*(u(2)/m);
end
