    
    clear all
    close all
    clc
    %% Data
m   = 1;          % massa [kg]
Cd  = 0.15;       % coefficiente di resistenza
t0  = 0;          % tempo iniziale [s]
tf  = 1.6;         % tempo finale [s]
 global nx nu N xc yc r h FigTraj

% Stati: [x; y; theta; v]
pos_i = [0; 0; 0];  
pos_f = [1;1; pi/3];

vel_i = 0;                 % velocità iniziale
vel_f = 0;                 % velocità finale

x_i = [pos_i; vel_i];
x_f = [pos_f; vel_f];

% max state vector used for normalization/scaling
xmax_vec = [pos_f(1:2); 2*pi/3;1.5];
% max control vector used for normalization/scaling
umax_vec = [1, 1]';
%% ----------------------- Pesi e vincoli del costo -----------------------
w1 = 0.2;
w2 = 0.3;
%final state weight
P = diag(1./xmax_vec.^2);
P = P.*diag([150 150 20 0]);
%optimal control weight
R = diag(1./umax_vec.^2);
R = R.*diag([w1 w2]);

Q = 0.05;                       % peso su penalizzazione velocità

alpha = 2;                    % peso vincolo morbido (soft constraint)
sigma = 0.01;                     % parametro di “softness”
xc = 0.65;
yc = 0.65;
r = 0.15;
            
   %state derivative 
  dx = @(x,u)[ x(4)*cos(x(3)); 
               x(4)*sin(x(3));
                 u(2);                                  
                (1/m)*(u(1) - Cd*(x(4))^2)];
    
    % Jacobian of the dynamics
    fx = @(x,u) [0 , 0 ,  -x(4)*sin(x(3))  ,  cos(x(3));
         0,  0,   x(4)*cos(x(3))  ,   sin(x(3))  ;
         0 ,  0  ,    0     ,            0  ;
         0 ,  0  ,  0 ,               -2*Cd*x(4)/m];
    
    
    fu = @(x,u) [0 , 0;
                0  ,  0;
                0  ,  1;
                1/m  ,  0];
    
       
    % Cost function
    Soft_cost_fun = @(x,u) alpha*exp((r^2 - (x(1)-xc).^2 - (x(2)-yc).^2)/sigma);
    
    L = @(x,u)  0.5*(u(1).*(R(1,1)*u(1)) + u(2).*(R(2,2)*u(2))) ...
           + Q*(x(4).^3) + Soft_cost_fun(x,u);    % Running cost
    
    p = @(x) 0.5*(x - x_f)'*P*(x - x_f) ; % Final cost
    
    a =  @(x,u) -alpha*(2*(x(1)-xc)/sigma)*exp((r^2 - (x(1)-xc).^2 -(x(2)-yc).^2)/sigma);
    b =  @(x,u) -alpha*(2*(x(2)-yc)/sigma)*exp((r^2 - (x(1)-xc)^2 -(x(2)-yc).^2)/sigma);
    
    
    Lx = @(x,u) [a(x,u) ,b(x,u) , 0, 3*Q*x(4).^2] ;
    
    Lu =@(x,u) [R(1,1)*u(1) , R(2,2)*u(2)];
    
    
    px =@(x) P*(x - x_f);
    
    %% Setup minimization problem
    
    nx = 4; % Number of states
    nu = 2; % number of inputs
    N = 301; % number of temporal steps
    h = tf/(N-1); % temporal discretization 
    
    
    options = optimoptions('fmincon', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'Display', 'iter', ...
    'MaxIterations', 40, ...
    'OutputFcn', @plotTrajectoryCallback, ...
    'PlotFcn', [], ...
    'OptimalityTolerance', 1e-5, ...    % 🔹 arresto se gradiente piccolo
    'StepTolerance', 1e-7, ...          % 🔹 arresto se il passo è piccolo
    'FunctionTolerance', 1e-5);         % 🔹 arresto se J varia poco

    %___ store the necessary stuff in a structure ____________________________%
    
    param.N = N;
    param.nu = nu;
    param.nx = nx;
    param.dx = dx;
    param.x_i = x_i;
    param.fx = fx;
    param.fu = fu;
    param.L = L;
    param.Lx = Lx;
    param.Lu = Lu;
    param.p = p;
    param.px = px;
    param.h = h;
    param.xc = xc;     % già definiti sopra
    param.yc = yc;
    param.r  = r;
    
    % Initial conditions for the minimization problem
    %z0 = ones(N*(nx + nu) + nx,1);
    % ... dopo aver impostato param, x_i, x_f, xc,yc,r
    z0 = build_initial_guess(param, x_i, x_f,m,Cd);
    
    %___ Define objective function and constraints ___________________________%
    
    ObjFun = @(z) cost_and_grad(z,param);
    NLcon = @(z) con_and_grad(z,param);
    
    % linear inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = []; 
    lb = []; % lower bound
    ub = []; % upper bound
    
   

    % Minimize ObjFun s.t. NLcon
    [z,fval] = fmincon(ObjFun,z0,A,b,Aeq,beq,lb,ub,NLcon,options);
    

 