% corrected_script.m
clear
close all
clc

% Data
m = 1;
Cd = 0.35;
t0 = 0;
tf = 15;

pos_i = [0; 0; pi/4];
pos_f = [10; 10; 0];
vel_i = 0;
vel_f = 1;

z_i = [pos_i; vel_i];
z_f = [pos_f; vel_f];

% weight for cost function
R = [2, 0; 0, 2];
Q = 0.05;
% weight for final state
P = diag([200,200,1,1]);
% soft boundary condition
alpha = 100;   % penalità meno rigida
sigma = 0.1;  % più ampio effetto dell’ostacolo
r = 0.1;
xc = 5;
yc = 5;

% Time interval (control parameterization grid)
Nsegment = 5000;                      % Number of control knots (can reduce)
Tu = linspace(t0, tf, Nsegment);    % control time vector

%% Iterative Procedure
options = odeset('RelTol', 1e-4, 'AbsTol',[1e-4 1e-4 1e-4 1e-4]);

Nmax = 5000;                        % Maximum number of iterations
u = repmat([1; 0], 1, Nsegment);    % initial control guess (2 x Nsegment)
step = 5e-6;                        % step size for gradient descent
eps = 1e-2;                         % exit tolerance on gradient norm

J = nan(1,Nmax);

for ii = 1:Nmax

    % 1) forward integrate state with current control
    [Tz,Z] = ode45(@(t,z) stateEq(t,z,u,Tu,m,Cd), [t0 tf], z_i, options);

    % Ensure column shapes
    Z = reshape(Z, length(Tz), 4);

    % 2) terminal adjoint
    Z_tf = Z(end,:).';
    lmb_tf = P*(Z_tf - z_f);

    % 3) backward integrate adjoint (note: AdjointEq expects access to Z & Tz)
    [Tlmb,lmb] = ode45(@(t,lmb) AdjointEq(t,lmb, Z,Tz, Q, sigma, alpha, xc, yc, r, Cd, m), [tf t0], lmb_tf, options);
    % lmb is on Tlmb (descending). Flip to ascending t to match Tz
    lmb = flipud(lmb);
    Tlmb = flipud(Tlmb);

    % Interpolate adjoint onto Tz (common grid)
    lmb1 = interp1(Tlmb, lmb(:,1), Tz);
    lmb2 = interp1(Tlmb, lmb(:,2), Tz);
    lmb3 = interp1(Tlmb, lmb(:,3), Tz);
    lmb4 = interp1(Tlmb, lmb(:,4), Tz);
    LMB = [lmb1 lmb2 lmb3 lmb4]';  % 4 x Nt

    % 4) compute dH/du on Tz grid
    dH = dHdu(u, Tu, LMB, Tz, m, R);   % returns Nt x 2

    H_Norm = norm(dH, 'fro');

    % 5) compute cost J(ii) accurately on Tz grid
    % Interpolate control to Tz
    u1_Tz = interp1(Tu, u(1,:), Tz)';
    u2_Tz = interp1(Tu, u(2,:), Tz)';

    z1 = Z(:,1);
    z2 = Z(:,2);
    z3 = Z(:,3);
    z4 = Z(:,4);

    soft_bound = @(Z1, Z2) alpha*exp((r^2 - (Z1-xc)^2 -(Z2-yc)^2)/sigma);
    J(ii) = 0.5*(Z_tf - z_f).'*P*(Z_tf - z_f);
    integrand = 0.5*(u1_Tz.^2 + u2_Tz.^2)*R(1,1) + Q*z4.^3 + soft_bound(z1,z2);
    J(ii) = J(ii) + trapz(Tz, integrand);


    % integrate with trapezoidal rule over Tz
    running_cost = trapz(Tz, integrand);
    final_cost = 0.5*(Z_tf - z_f).' * P * (Z_tf - z_f);
    J(ii) = final_cost + running_cost;

    % show progress every some iterations
    if mod(ii,50)==1
        fprintf('Iter %d: J = %.4f, ||dH|| = %.6f\n', ii, J(ii), H_Norm);
    end

    % 6) check convergence
    if H_Norm < eps
        disp(['Converged at iteration ', num2str(ii), '. Final cost: ', num2str(J(ii))]);
        break;
    end

    % 7) update control: need dH on Tz mapped back to Tu (control grid)
    % Interpolate dH (Nt x 2) to Tu
    dH1_onTu = interp1(Tz, dH(1,:), Tu);
    dH2_onTu = interp1(Tz, dH(2,:), Tu);
    dH_onTu = [dH1_onTu; dH2_onTu];

    % Gradient step (element-wise)
    u = u - step * dH_onTu;

    % Optional: clamp controls if physically bounded (example)
    % u(1,:) = max(min(u(1,:),u1_max), u1_min);
    % u(2,:) = max(min(u(2,:),u2_max), u2_min);

end

% final outputs: Tz, Z, u, J
disp('Done.');

%% Subfunctions (same names but fixed/consistent vectorization)

function u_new = IterControl(dH,tx,u,tu,step)
    % NOT used in corrected loop; left for compatibility
    dH1 = interp1(tx,dH(:, 1),tu);
    dH2 = interp1(tx,dH(:, 2),tu);
    dH_int = [dH1; dH2];
    u_new = u - step*dH_int;
end

function dH = dHdu(u, Tu, LMB, tz, m, R)
    % Interpolate control to tz (row arrays)
    u1 = interp1(Tu, u(1, :), tz);
    u2 = interp1(Tu, u(2, :), tz);

    % L_u: for each time sample, cost derivative R*u
    % Compute R*[u1; u2] for all times: (2 x Nt) result
    Umat = [u1'; u2'];           % 2 x Nt
    RU = R * Umat;            % 2 x Nt
    L_u = RU.';               % Nt x 2

    % LMB is 4 x Nt ; LMB' (Nt x 4) * B (4 x 2) => Nt x 2
    B  = [ 0,  0;
           0,  0;
           0,  1;
          1/m, 0];
    dH2 = L_u + (LMB.' * B);   % Nt x 2
    dH = (R*[u1'; u2']) + B'*LMB;

end

function dz = stateEq(t,z,u,Tu,m,Cd)
    dz = zeros(4,1);
    % interpolate control at time t
    u1 = interp1(Tu, u(1, :), t);
    u2 = interp1(Tu, u(2, :), t);

    dz(1) = z(4)*cos(z(3));
    dz(2) = z(4)*sin(z(3));
    dz(3) = u2;
    dz(4) = 1/m*(u1 - Cd*z(4)^2);
end

function dlmb = AdjointEq(t,lmb, Z,Tz, Q, sigma, alpha, xc, yc, r, Cd, m)
    % Interpolate states at time t
    z4 = interp1(Tz, Z(:,4), t);
    z3 = interp1(Tz, Z(:,3), t);
    z2 = interp1(Tz, Z(:,2), t);
    z1 = interp1(Tz, Z(:,1), t);

    A = [0   0  -z4*sin(z3)  cos(z3);
         0   0   z4*cos(z3)  sin(z3);
         0   0     0           0;
         0   0     0        -2*Cd*z4/m];

    % derivatives of soft bound wrt x and y
    expo = exp((r^2 - (z1-xc)^2 -(z2-yc)^2)/sigma);
    a = -alpha*(2*(z1-xc)/sigma)*expo;
    b = -alpha*(2*(z2-yc)/sigma)*expo;

    L_x = [a b 0 3*Q*z4^2];

    % dlmb = -(L_x + lmb' * A)'  (same as your original)
    dlmb = -(L_x + (lmb' * A))';
end
