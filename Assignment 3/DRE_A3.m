function [OutVect] = DRE(t,p,Q,R,Z,Tz,u,Tu)

% evaluation of the state at time @t
Nstates = size(Z,2);
zk = zeros(1,Nstates);
for ii = 1:Nstates
    zk(:,ii) = interp1(Tz,Z(:,ii),t);
end

% evaluation of the control action at time @t
Nu = size(u,2);
uk = zeros(1,Nu);
for ii = 1:Nu
    uk(:,ii) = interp1(Tu,u(:,ii),t);
end
m = 1;
Cd = 0.15;

A = [0 0 -zk(4)*sin(zk(3))  cos(zk(3));
         0 0  zk(4)*cos(zk(3))  sin(zk(3));
         0 0  0           0;
         0 0  0       -2*Cd*zk(4)/m];

B = [0 0; 0 0; 0 1; 1/m 0];

% transformation vector -> matrix
P = zeros(Nstates,Nstates);
P(1:end) = p(1:end);

% DRE
Out = -(A.'*P + P*A - P*B*R^-1*B.'*P + Q);

% transformation matrix -> vector
OutVect = Out(1:end)';