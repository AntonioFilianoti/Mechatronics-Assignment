function [cost, grad] = cost_and_grad(z, param)
% =========================================================================
% COST_AND_GRAD
% -------------------------------------------------------------------------
% Scopo:
%   Calcola la funzione obiettivo totale e il suo gradiente rispetto al
%   vettore decisionale z, per il metodo diretto di controllo ottimo.
%
% Descrizione:
%   Il costo è definito come:
%       J = Σ [ h * L(x_k, u_k) ] + p(x_{N+1})
%   dove L è il costo "istantaneo" (running cost),
%   e p è il costo finale (terminal cost).
%
% Input:
%   z      - vettore decisionale (stati + controlli concatenati)
%   param  - struttura parametri (L, Lx, Lu, p, px, N, nx, nu, h)
%
% Output:
%   cost - valore scalare della funzione obiettivo
%   grad - gradiente colonna (stessa dimensione di z)
% =========================================================================

%% ---------------------- Estrazione parametri ----------------------------
N  = param.N;
nx = param.nx;
nu = param.nu;
h  = param.h;

L  = param.L;   % costo istantaneo
p  = param.p;   % costo finale

%% --------------------- Estrazione stati e controlli ---------------------
x = zeros(nx, N+1);
u = zeros(nu, N);

for k = 1:N+1
    idx = (k-1)*(nx+nu) + (1:nx);
    x(:,k) = z(idx);
end

for k = 1:N
    idx = (k-1)*(nx+nu) + nx + (1:nu);
    u(:,k) = z(idx);
end

%% ------------------------ Calcolo del costo -----------------------------
cost = 0;
for k = 1:N
    cost = cost + h * L(x(:,k), u(:,k));  % costo istantaneo integrato
end
cost = cost + p(x(:,end));                % costo finale

%% ---------------------- Calcolo del gradiente ---------------------------
if nargout > 1
    Lx = param.Lx;
    Lu = param.Lu;
    px = param.px;

    grad = zeros(size(z));   % gradiente complessivo

    % contributo di L(x,u)
    for k = 1:N
        idx_x = (k-1)*(nx+nu) + (1:nx);
        idx_u = (k-1)*(nx+nu) + nx + (1:nu);

        grad(idx_x) = grad(idx_x) + h * Lx(x(:,k), u(:,k))';
        grad(idx_u) = grad(idx_u) + h * Lu(x(:,k), u(:,k))';
    end

    % contributo del costo finale p(x_N)
    grad(end - nx + 1:end) = grad(end - nx + 1:end) + px(x(:,end));
end

end
