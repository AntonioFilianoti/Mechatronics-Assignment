function [dlmb] = AdjointEq(t,lmb, Z,Tz, Q, sigma, alpha, xc, yc, r, Cd, m)
    dlmb = zeros(4,1);
    z4 = interp1(Tz,Z(:,4),t);   % Interploate the state varialbes
    z3 = interp1(Tz,Z(:,3),t);   % Interploate the state varialbes
    z2 = interp1(Tz,Z(:,2),t);   % Interploate the state varialbes
    z1 = interp1(Tz,Z(:,1),t);   % Interploate the state varialbes

    A = [0   0  -z4*sin(z3)  cos(z3);
         0   0   z4*cos(z3) sin(z3);
         0   0     0           0;
         0   0     0        -2*Cd*z4/m];
    

    expo = (r^2 - (z1-xc)^2 -(z2-yc)^2)/sigma;
    expo = min(expo, 1000);  % previene overflow numerico
    exp_term = exp(expo);

    a = -alpha*(2*(z1-xc)/sigma)*exp_term;
    b = -alpha*(2*(z2-yc)/sigma)*exp_term;

    L_x = [a b 0 3*Q*z4^2];

    % lambda equation
    dlmb = -(L_x + lmb' * A)';
end
