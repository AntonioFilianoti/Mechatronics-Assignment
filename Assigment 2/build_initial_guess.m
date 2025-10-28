function z0 = build_initial_guess(param, x_i, x_f, m, Cd)
% =========================================================================
% BUILD_INITIAL_GUESS
% -------------------------------------------------------------------------
% Scopo:
%   Costruisce una stima iniziale (z0) coerente con la dinamica del sistema
%   da utilizzare come punto di partenza per l’ottimizzazione con fmincon.
%
% Descrizione:
%   Genera una traiettoria iniziale “ragionevole” che collega lo stato
%   iniziale x_i e quello finale x_f, approssimando posizione, velocità e
%   controlli in modo coerente con le equazioni dinamiche del veicolo.
%
% Input:
%   param  - struttura con i parametri del problema (N, nx, nu, h, dx, ecc.)
%   x_i    - stato iniziale [x; y; θ; v]
%   x_f    - stato finale desiderato [x; y; θ; v]
%   m, Cd  - massa e coefficiente di resistenza
%
% Output:
%   z0     - vettore colonna contenente tutti gli stati e i controlli
%             concatenati (lunghezza: N*(nx+nu) + nx)
% =========================================================================

%% ------------------------- Estrazione parametri -------------------------
N  = param.N;   % numero di passi temporali
nx = param.nx;  % dimensione stato
nu = param.nu;  % dimensione controllo
h  = param.h;   % passo temporale (Δt)
dx = param.dx;  % funzione dinamica: x_dot = f(x,u)

%% ----------------------- Preallocazione variabili -----------------------
% x: matrice degli stati [nx x (N+1)] → include lo stato finale
% u: matrice dei controlli [nu x N]
x = zeros(nx, N+1);
u = zeros(nu, N);

%% ---------------- Interpolazione posizioni e velocità -------------------
% Interpola linearmente le coordinate x e y da iniziale a finale
for i = 1:2
    x(i, :) = linspace(x_i(i), x_f(i), N+1);
end

% Calcolo dell’orientamento (theta) come direzione tangente alla traiettoria
dx_path = diff(x(1, :));   % differenza lungo x
dy_path = diff(x(2, :));   % differenza lungo y
phi = atan2(dy_path, dx_path);  % angolo di direzione istantanea
phi = [phi, phi(end)];           % estendi per avere N+1 punti

% Interpolazione della velocità scalare da v_i a v_f
v = linspace(x_i(4), x_f(4), N+1);

% Inserisci phi e v nel vettore di stato
x(3, :) = phi;   % orientamento
x(4, :) = v;     % velocità

%% ------------------- Calcolo dei controlli coerenti ---------------------
% Si ricavano i controlli “coerenti” con la traiettoria interpolata.
% L’obiettivo è stimare u in modo che la dinamica non venga violata troppo.

for k = 1:N
    % Controllo angolare u2: variazione dell’orientamento
    % u2 = dφ/dt ≈ (φ(k+1) - φ(k)) / h
    u(2, k) = (phi(k+1) - phi(k)) / h;

    % Controllo di trazione u1: dalla dinamica di v_dot
    % v_dot = (v(k+1) - v(k)) / h  →  F = m*v_dot + Cd*v^2
    vdot = (v(k+1) - v(k)) / h;
    u(1, k) = m * vdot + Cd * v(k)^2;
end

%% ---------------- Correzione posizione con integrazione -----------------
% Integra la dinamica “vera” del sistema per correggere le posizioni.
% In questo modo si ottiene una traiettoria iniziale dinamicamente più coerente.
for k = 1:N
    x(:, k+1) = x(:, k) + h * dx(x(:, k), u(:, k));  % integrazione di Eulero
end

%% ---------------- Riallineamento alla posizione finale ------------------
% Dopo l’integrazione, la posizione finale potrebbe non coincidere esattamente
% con x_f. Si effettua quindi una correzione di scala lineare.

% Fattore di scala per correggere la traiettoria
scale = (x_f(1) - x(1, 1)) / (x(1, end) - x(1, 1) + eps);

% Riscala le coordinate x e y per forzare l’allineamento alla destinazione
x(1, :) = x(1, 1) + scale * (x(1, :) - x(1, 1));
x(2, :) = x(2, 1) + scale * (x(2, :) - x(2, 1));

% Allinea orientamento e velocità finali ai valori desiderati
x(3, end) = x_f(3);
x(4, end) = x_f(4);

%% -------------------- Impacchettamento nel vettore z0 -------------------
% z0 è il vettore che contiene in sequenza:
% [x(1); u(1); x(2); u(2); ... ; x(N); u(N); x(N+1)]

z0 = zeros(N * (nx + nu) + nx, 1);  % vettore finale

for k = 1:N
    base = (k - 1) * (nx + nu);
    z0(base + (1:nx))      = x(:, k);   % stato al passo k
    z0(base + nx + (1:nu)) = u(:, k);   % controllo al passo k
end

% Aggiungi l’ultimo stato (x_{N+1})
z0(N * (nx + nu) + (1:nx)) = x(:, end);

end
