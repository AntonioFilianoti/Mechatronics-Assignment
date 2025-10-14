function [dz] = stateEq(t,z,u,Tu,m,Cd)
    
    u1 = interp1(Tu,u(1, :),t);
    u2 = interp1(Tu,u(2, :),t);
   
    
    dz = zeros(4, 1);

    dz(1) = z(4)*cos(z(3));
    dz(2) = z(4)*sin(z(3));
    dz(3) = u2;
    dz(4) = 1/m*(u1-Cd*z(4)^2);
end