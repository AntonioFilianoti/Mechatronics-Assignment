function u = IterControl_new(dH, Tz, u_old, Tu, step)
    dH1 = interp1(Tz,dH(1,:), Tu);
    dH2 = interp1(Tz,dH(2,:), Tu);

    dH = [dH1; dH2];

    u = u_old - step*dH;
end