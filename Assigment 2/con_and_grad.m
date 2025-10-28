function [c, ceq, gradc, gradceq] = con_and_grad(z, param)
% =========================================================================
% CON_AND_GRAD
% -------------------------------------------------------------------------
% Scopo:
%   Definisce i vincoli di uguaglianza (dinamica e condizione iniziale)
%   e disuguaglianza (evitamento ostacolo) per il problema di ottimizzazione
%   non lineare diretto (NLP), insieme ai rispettivi gradienti analitici.
%
% Descrizione:
%   I vincoli di uguaglianza assicurano che la dinamica del sistema sia
%   rispettata fra stati consecutivi (integrazione di Eulero esplicito).
%   I vincoli di disuguaglianza modellano un ostacolo circolare nel piano.
%
% Input:
%   z      - vettore decisionale contenente stati e controlli concatenati
%   param  - struttura parametri (N, nx, nu, dx, fx, fu, x_i, xc, yc, r)
%
% Output:
%   c        - vincoli di disuguaglianza (c <= 0)
%   ceq      - vincoli di uguaglianza (ceq = 0)
%   gradc    - Jacobiano di c rispetto a z
%   gradceq  - Jacobiano di ceq rispetto a z
% =========================================================================

%% ---------------------- Estrazione parametri ----------------------------
N   = param.N;     % numero di intervalli temporali
nx  = param.nx;    % dimensione stato (4)
nu  = param.nu;    % dimensione controllo (2)
h   = param.h;     % passo temporale
dx  = param.dx;    % funzione dinamica f(x,u)
fx  = param.fx;    % Jacobiano ∂f/∂x
fu  = param.fu;    % Jacobiano ∂f/∂u
x_i = param.x_i;   % stato iniziale noto

% Parametri ostacolo
xc  = param.xc;
yc  = param.yc;
r   = param.r;

%% --------------------- Estrazione di x e u da z ------------------------
x = zeros(nx, N+1);
u = zeros(nu, N);

for k = 1:N+1
    base = (k-1)*(nx+nu);
    x(:,k) = z(base + (1:nx));   % estrai stato
end

for k = 1:N
    base = (k-1)*(nx+nu);
    u(:,k) = z(base + nx + (1:nu)); % estrai controllo
end

%% ==================== VINCOLI DI UGUAGLIANZA ============================
% Ceq include:
% (1) vincoli di dinamica  → x_{k+1} = x_k + h*f(x_k,u_k)
% (2) vincolo iniziale     → x_1 = x_i

ceq = zeros(nx*N + nx, 1);         % vettore di vincoli
nz = numel(z);                     % dimensione vettore z
gradceq = sparse(nz, nx*N + nx);   % Jacobiano dei vincoli di uguaglianza

I = eye(nx);  % matrice identità utile per le derivate

% --- (1) Dinamica discretizzata ---
for k = 1:N
    idx_c = (k-1)*nx + (1:nx);          % posizione dei vincoli
    fk = dx(x(:,k), u(:,k));            % f(x_k,u_k)
    
    % vincolo: x_{k+1} - x_k - h*f(x_k,u_k) = 0
    ceq(idx_c) = x(:,k+1) - x(:,k) - h*fk;
    
    % derivate parziali del vincolo rispetto a z
    Ak = -I - h*fx(x(:,k), u(:,k));   % ∂/∂x_k
    Bk = -h*fu(x(:,k), u(:,k));       % ∂/∂u_k
    Ck =  I;                          % ∂/∂x_{k+1}
    
    zxk   = (k-1)*(nx+nu) + (1:nx);
    zuk   = (k-1)*(nx+nu) + nx + (1:nu);
    zxkp1 = (k)*(nx+nu) + (1:nx);
    
    % Inserimento trasposto per coerenza con fmincon
    gradceq(zxk,   idx_c) = Ak';
    gradceq(zuk,   idx_c) = Bk';
    gradceq(zxkp1, idx_c) = Ck';
end

% --- (2) Condizione iniziale ---
idx_c0 = nx*N + (1:nx);
ceq(idx_c0) = x(:,1) - x_i;  % x_1 - x_i = 0

% Derivata rispetto a x_1 (le altre variabili non compaiono)
zx1 = 0*(nx+nu) + (1:nx);
gradceq(zx1, idx_c0) = I;

%% ==================== VINCOLI DI DISUGUAGLIANZA =======================
% Evitamento ostacolo:
%   c_k = r^2 - (x_k - xc)^2 - (y_k - yc)^2  <= 0
%   se c_k > 0 → violazione: punto all’interno dell’ostacolo

c = zeros(N+1, 1);
gradc = sparse(nz, N+1);

for k = 1:N+1
    Xk = x(:,k);
    dxk = Xk(1) - xc;    % differenza su x
    dyk = Xk(2) - yc;    % differenza su y

    % vincolo
    c(k) = r^2 - dxk^2 - dyk^2;

    % gradiente (solo x,y influenzano il vincolo)
    g = zeros(nx,1);
    g(1) = -2*dxk;
    g(2) = -2*dyk;

    zxk = (k-1)*(nx+nu) + (1:nx);
    gradc(zxk, k) = g;
end

end
