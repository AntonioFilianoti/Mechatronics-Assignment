clear
close all
clc

% Data
m = 1;
Cd = 0.35;
t0 = 0;
tf = 10;

pos_i = [0; 0; pi/4];
pos_f = [10; 10; 0];
vel_i = 0;
vel_f = 1;

z_i = [pos_i; vel_i];
z_f = [pos_f; vel_f];


% weight for cost function
R = [20, 0; 0, 20];
Q = 0.05;
% weight for final state
P = diag([100,100,1,1]);
% soft boundary condition
alpha = 10000;
sigma = 10;
r = 1;
xc = 5;
yc = 5;

theta = linspace(0, 2*pi, 200); % più punti per un cerchio più liscio
x_cons = xc + r .* cos(theta);
y_cons = yc + r .* sin(theta);

figure('Position', [100 100 900 700]); % <-- figura più grande (x,y,width,height)
plot(x_cons, y_cons, 'LineWidth', 1.8);
grid on; grid minor;
hold on;
scatter(pos_i(1), pos_i(2), 60, 'filled');
scatter(pos_f(1), pos_f(2), 60, 'filled');

axis equal;

% Limiti più larghi (ad esempio, un margine del 20% rispetto al raggio)
xlim([-10, 20]);
ylim([-10, 20]);

xlabel('x');
ylabel('y');
legend('Constraint', 'Start', 'Finish');



% Time interval
Nsegment = 1000;                      % Number of time intervals considered
Tu = linspace(t0, tf, Nsegment);    % discretize time

%% Iterative Procedure
options = odeset('RelTol', 1e-4, 'AbsTol',[1e-4 1e-4 1e-4 1e-4]);

Nmax = 5000;                        % Maximum number of iterations
u = [1; 0]*ones(1,Nsegment);               % guessed initial control  u = 1
step = 0.00005;                       % speed of control adjustment                      % speed of control adjustment
eps = 1e-2;                         % Exit tollerance condition

for ii = 1:Nmax

    % 1) start with assumed control u and move forward
    [Tz,Z] = ode45(@(t,z) stateEq(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);
    % [Tz,Z] = ode45(@(t,z) stateEq(t,z,u,Tu,m,Cd), Tu, z_i, options);

    % 2) Move backward to get the adjoint vectos trajectory
    Z_tf = Z(end,:).';
    lmb_tf = P*(Z_tf - z_f);

    [Tlmb,lmb] = ode45(@(t,lmb) AdjointEq(t,lmb, Z,Tz, Q, sigma, alpha, xc, yc, r, Cd, m), [t0 tf], lmb_tf, options);
    lmb = flipud(lmb);
    Tlmb = flipud(Tlmb);
    lmb1 = lmb(:,1);
    lmb1 = interp1(Tlmb,lmb1,Tz);
    lmb2 = lmb(:,2);
    lmb2 = interp1(Tlmb,lmb2,Tz);
    lmb3 = lmb(:,3);
    lmb3 = interp1(Tlmb,lmb3,Tz);
    lmb4 = lmb(:,4);
    lmb4 = interp1(Tlmb,lmb4,Tz);

    LMB = [lmb1 lmb2 lmb3 lmb4]';
LMB = interp1(Tz, LMB', Tu, 'linear', 'extrap')';

    dH = dHdu(u, Tu, LMB, Tz, m, R);

    H_Norm = norm(dH, 'fro');
    if H_Norm >= 1e10
        H_Norm = 1e9;
        dH = dH./1e6;
    end

    Z4 = interp1(Tz,Z(:,4),Tu);   % Interploate the state varialbes
    Z3 = interp1(Tz,Z(:,3),Tu);   % Interploate the state varialbes
    Z2 = interp1(Tz,Z(:,2),Tu);   % Interploate the state varialbes
    Z1 = interp1(Tz,Z(:,1),Tu);


    soft_bound = @(Z1, Z2) alpha*exp((r^2 - (Z1-xc)^2 -(Z2-yc)^2)/sigma);
    % J(ii) = 0.5*(Z_tf - z_f).'*P*(Z_tf - z_f);
    % dt = Tz(2) - Tz(1);  % passo costante
    % for k = 1:length(Tz)
    %     J(ii) = J(ii) + dt*(0.5*(u(:,k)'*R*u(:,k)) + Q*Z4(k)^3 + soft_bound(Z1(k),Z2(k)));
    % end

end

% disp(['Final cost: ',num2str(J(ii)),' [-]'])
% disp(' ')
plot(Z1,Z2)