function z0 = build_initial_guess(param, x_i, x_f,m,Cd)
%BUILD_INITIAL_GUESS_DYNAMIC crea una z0 coerente con la dinamica
% per qualsiasi condizione iniziale e finale (senza parametri manuali)

% === Estrai parametri
N  = param.N;
nx = param.nx;
nu = param.nu;
h  = param.h;
dx = param.dx;


% === Preallocazione
x = zeros(nx, N+1);
u = zeros(nu, N);

% === Interpola posizioni e velocità
for i = 1:2
    x(i,:) = linspace(x_i(i), x_f(i), N+1);
end

% orientamento φ: direzione tangente alla traiettoria
dx_path = diff(x(1,:));
dy_path = diff(x(2,:));
phi = atan2(dy_path, dx_path);
phi = [phi, phi(end)]; % per lunghezza N+1

% velocità scalare (lineare da v_i a v_f)
v = linspace(x_i(4), x_f(4), N+1);

% salva nel vettore di stato
x(3,:) = phi;
x(4,:) = v;

% === Controlli coerenti
for k = 1:N
    % u2 = dφ/dt
    u(2,k) = (phi(k+1) - phi(k)) / h;

    % u1 = m * vdot + Cd * v^2
    vdot = (v(k+1) - v(k)) / h;
    u(1,k) = m*vdot + Cd*v(k)^2;
end

% === Correggi le posizioni integrando la dinamica
for k = 1:N
    x(:,k+1) = x(:,k) + h * dx(x(:,k), u(:,k));
end

% === riallinea il finale alla posizione richiesta
scale = (x_f(1) - x(1,1)) / (x(1,end) - x(1,1) + eps);
x(1,:) = x(1,1) + scale * (x(1,:) - x(1,1));
x(2,:) = x(2,1) + scale * (x(2,:) - x(2,1));
x(3,end) = x_f(3);  % allinea phi finale
x(4,end) = x_f(4);  % allinea v finale

% === impacchetta
z0 = zeros(N*(nx+nu)+nx,1);
for k = 1:N
    base = (k-1)*(nx+nu);
    z0(base + (1:nx))      = x(:,k);
    z0(base + nx + (1:nu)) = u(:,k);
end
z0(N*(nx+nu) + (1:nx)) = x(:,end);
