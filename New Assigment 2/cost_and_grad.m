function [cost, grad] = cost_and_grad(z, param)
% COST_AND_GRAD  Objective and gradient for NLP direct method
%
%   [cost, grad] = cost_and_grad(z, param)
%
%   Computes the total cost:
%       J = sum_k h*L(x_k, u_k) + p(x_N)
%   and its gradient w.r.t. decision vector z
%
%   where z = [x_1; u_1; x_2; u_2; ... ; x_N; u_N; x_{N+1}]
%
%   Inputs:
%       param.L, param.Lx, param.Lu, param.p, param.px
%       param.N, param.nx, param.nu, param.h
%
%   Output:
%       cost - scalar objective value
%       grad - gradient column vector (same length as z)

% === Unpack parameters
N  = param.N;
nx = param.nx;
nu = param.nu;
h  = param.h;
L  = param.L;
p  = param.p;

% === Extract states and controls from z
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

% === Cost function
cost = 0;
for k = 1:N
    cost = cost + h * L(x(:,k), u(:,k));
end
cost = cost + p(x(:,end));

% === Gradient of the cost function
if nargout > 1
    Lx = param.Lx;
    Lu = param.Lu;
    px = param.px;

    grad = zeros(size(z));

    for k = 1:N
        % indices in z for x_k and u_k
        idx_x = (k-1)*(nx+nu) + (1:nx);
        idx_u = (k-1)*(nx+nu) + nx + (1:nu);

        grad(idx_x) = grad(idx_x) + h * Lx(x(:,k), u(:,k))';
        grad(idx_u) = grad(idx_u) + h * Lu(x(:,k), u(:,k))';
    end

    % terminal contribution
    grad(end - nx + 1:end) = grad(end - nx + 1:end) + px(x(:,end));
end
end
