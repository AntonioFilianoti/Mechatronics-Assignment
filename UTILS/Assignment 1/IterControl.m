function u_new = IterControl(dH,tx,u,tu,step)
% interploate dH/du
dH1 = interp1(tx,dH(1, :),tu);
dH2 = interp1(tx,dH(2, :),tu);
dH_int = [dH1; dH2];
u_new = u - step*dH_int;
end