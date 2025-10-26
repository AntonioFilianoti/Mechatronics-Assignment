function [dH ] = dHdu_new(u, Tu,lmb4, lmb3,Tz, m, R)
    u1 = interp1(Tu,u(1,:),Tz);
    u2 = interp1(Tu,u(2,:),Tz);
    u = [u1 u2]';


    %Lu = R*u;

    B = [0 0; 0 0; 0 1; 1/m 0];

   % dH_1 = Lu + B'*LMB;
    dH = [u1*R(1,1) + (lmb4/m), lmb3 + u2*R(2,2)]';
end