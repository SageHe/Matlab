function[dydt] = hw3ode(t,y)
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4
%define perturbations to the system;
d_p = y(1);
d_q = y(2);
d_r = y(3);
d_phi = y(4);
d_theta = y(5);
d_u = y(6);
d_v = y(7);
d_w = y(8);

d_Lc = 0;
d_Mc = 0;
d_Nc = 0;
d_Zc = 0;

deltaq_dot = (1/I_y)*d_Mc;
deltap_dot = (1/I_x)*d_Lc;
deltar_dot = (1/I_z)*d_Nc;
deltatheta_dot = d_q;
deltaphi_dot = d_p;
deltaU_dot = -g*d_theta;
deltaV_dot = g*d_phi;
deltaW_dot = (1/m)*d_Zc;

dydt(1) = deltaq_dot;
dydt(2) = deltap_dot;
dydt(3) = deltar_dot;
dydt(4) = deltaphi_dot;
dydt(5) = deltatheta_dot;
dydt(6) = deltaU_dot;
dydt(7) = deltaV_dot;
dydt(8) = deltaW_dot;

dydt = dydt';

