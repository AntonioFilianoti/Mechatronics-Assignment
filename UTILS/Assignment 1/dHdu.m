function dH = dHdu(u, Tu, LMB, tz, m, R)
% interploate the control
% u1 = interp1(Tu,u(1, :),tz)';
% u2 = interp1(Tu,u(2, :),tz)';
u1 = u(1,:);
u2 = u(2,:);


L_u = (R*[u1; u2]);

B  = [ 0,  0;
       0,  0;
       0,  1;
      1/m, 0];

dH = L_u + B'*LMB;
end