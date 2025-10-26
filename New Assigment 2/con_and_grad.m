function [c, ceq, gradc, gradceq] = con_and_grad(z, param)
% NLP direct method - equality (dynamics + initial) and inequality (obstacle) constraints
% with analytic Jacobians for fmincon.

% -------- unpack
N   = param.N;    % number of grid intervals (we have N+1 states and N controls)
nx  = param.nx;   % 4  -> [x;y;phi;v]
nu  = param.nu;   % 2  -> [u1;u2]
h   = param.h;
dx  = param.dx;   % f(x,u)
fx  = param.fx;   % df/dx (nx x nx)
fu  = param.fu;   % df/du (nx x nu)
x_i = param.x_i;  % initial state

% Obstacle parameters (center and radius)
xc  = param.xc;   
yc  = param.yc;
r   = param.r;

% -------- extract x,u from z 
x = zeros(nx, N+1);  u = zeros(nu, N);
for k = 1:N+1
    base = (k-1)*(nx+nu);
    x(:,k) = z(base + (1:nx));
end
for k = 1:N
    base = (k-1)*(nx+nu);
    u(:,k) = z(base + nx + (1:nu));
end

% =======================
% Equality constraints
% =======================
% (1) Dynamics:  x_{k+1} - x_k - h f(x_k,u_k) = 0   for k=1..N
% (2) Initial:   x_1 - x_i = 0
ceq = zeros(nx*N + nx, 1);

% Helper to place blocks into the big gradient matrix
nz = numel(z);
gradceq = sparse(nz, nx*N + nx);

I = eye(nx);

% (1) dynamics constraints
for k = 1:N
    idx_c = (k-1)*nx + (1:nx);               % vincoli per step k
    fk    = dx(x(:,k), u(:,k));              % f(x_k,u_k)

    % vincolo di dinamica
    ceq(idx_c) = x(:,k+1) - x(:,k) - h*fk;

    % Jacobiani
    Ak = -I - h*fx(x(:,k), u(:,k));          % ∂/∂x_k
    Bk =     - h*fu(x(:,k), u(:,k));         % ∂/∂u_k
    Ck =  I;                                 % ∂/∂x_{k+1}

    zxk   = (k-1)*(nx+nu) + (1:nx);
    zuk   = (k-1)*(nx+nu) + nx + (1:nu);
    zxkp1 = (k  )*(nx+nu) + (1:nx);

   % assegnazione trasposta ---
    gradceq(zxk , idx_c)   = Ak';
    gradceq(zuk , idx_c)   = Bk';
    gradceq(zxkp1, idx_c)  = Ck';
end


% (2) initial condition x_1 = x_i
idx_c0 = nx*N + (1:nx);
ceq(idx_c0) = x(:,1) - x_i;

zx1 = 0*(nx+nu) + (1:nx);
gradceq(zx1, idx_c0) = I;

% =======================
% Inequality constraints (obstacle)
% c <= 0
% c_k = r^2 - (x_k - xc)^2 - (y_k - yc)^2   for k=1..N+1
% =======================
c = zeros(N+1, 1);
gradc = sparse(nz, N+1);

for k = 1:N+1
    Xk = x(:,k);
    dxk = Xk(1) - xc;   % x - xc
    dyk = Xk(2) - yc;   % y - yc

    c(k) = r^2 - dxk^2 - dyk^2;

    % gradient wrt z: only depends on x_k entries (x,y)
    g = zeros(nx,1);
    g(1) = -2*dxk;
    g(2) = -2*dyk;
    % phi and v have zero derivatives for the obstacle
    zxk = (k-1)*(nx+nu) + (1:nx);  % valid also for k=N+1
    gradc(zxk, k) = g;
end

end
