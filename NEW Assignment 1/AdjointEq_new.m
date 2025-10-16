function dlmb = AdjointEq_new(t,lmb,Z,Tz,Q,sigma,alpha,xc,yc,r,Cd,m, lmb_tf )

    dlmb = zeros(4,1);
    z1 = interp1(Tz,Z(:,1), t);
    z2 = interp1(Tz,Z(:,2), t);
    z3 = interp1(Tz,Z(:,3), t);
    z4 = interp1(Tz,Z(:,4), t);

    % u1 = interp1(Tu,u(1,:),t);
    % u2 = interp1(Tu,u(2,:),t);
    
   A = [0 0 -z4*sin(z3) cos(z3);
         0 0  z4*cos(z3) sin(z3);
         0 0  0          0;
         0 0  0       -2*Cd*z4/m];
 
    % Gradiente del vincolo morbido
    expo = (r^2 - (z1-xc)^2 - (z2-yc)^2)/sigma;
    exp_term = alpha*exp(expo);
    a = -(2*(z1-xc)/sigma)*exp_term;
    b = -(2*(z2-yc)/sigma)*exp_term;

    % Gradiente del Lagrangiano rispetto a z
    Lz = [a b 0 3*Q*z4^2];
    dlmb = -(Lz + (lmb'*A))';
    % dlmb = -Lz' - A.' * lmb;
end