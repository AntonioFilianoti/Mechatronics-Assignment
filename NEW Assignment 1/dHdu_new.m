function dH = dHdu_new(u, Tu, LMB, Tz, m, R)
    u1 = interp1(Tu,u(1,:),Tz);
    u2 = interp1(Tu,u(2,:),Tz);
    u = [u1 u2]';


    Lu = R*u;

    B = [0 0; 0 0; 0 1; 1/m 0];

    dH = Lu + B'*LMB;
end