function [dydt] = hw12ode(t,y)
global Alat_aug Blat_aug K
Alat_aug_cl = (Alat_aug + Blat_aug*K);
ydot = Alat_aug_cl*y;

dydt = ydot;
% dydt = dydt';
end