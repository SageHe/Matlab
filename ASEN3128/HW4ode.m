function[dydt] = hw3ode(t,y)
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 K1_roll K2_roll K1_pitch K2_pitch K_y
%define perturbations to the system;
d_p = y(1);
d_q = y(2);
d_r = y(3);
d_phi = y(4);
d_theta = y(5);
% d_psi = y(6);
d_u = y(6);
d_v = y(7);
d_w = y(8);
% d_x = y(10);
% d_y = y(11);
% d_z = y(12);

d_Lc = (-K1_roll*d_p - K2_roll*d_phi);
d_Mc = (-K1_pitch*d_q - K2_pitch*d_theta);
d_Nc = -K_y*d_r;
d_Zc = 0;

deltaq_dot = (1/I_y)*d_Mc;
deltap_dot = (1/I_x)*d_Lc;
deltar_dot = (1/I_z)*d_Nc;
deltatheta_dot = d_q;
deltaphi_dot = d_p;
deltaU_dot = -g*d_theta;
deltaV_dot = g*d_phi;
deltaW_dot = (1/m)*d_Zc;

dydt(1) = deltap_dot;
dydt(2) = deltaq_dot;
dydt(3) = deltar_dot;
dydt(4) = deltaphi_dot;
dydt(5) = deltatheta_dot;
% dydt(6) = deltapsi_dot;
dydt(6) = deltaU_dot;
dydt(7) = deltaV_dot;
dydt(8) = deltaW_dot;
% dydt(10) = d_u*cos(d_theta)*cos(d_psi) + d_v*(sin(d_phi)*sin(d_theta)*cos(d_psi) - cos(d_phi)*sin(d_psi)) + d_w*(cos(d_phi)*sin(d_theta)*cos(d_psi) + sin(d_phi)*sin(d_psi));
% dydt(11) = d_u*cos(d_theta)*sin(d_psi) + d_v*(sin(d_phi)*sin(d_theta)*sin(d_psi) + cos(d_phi)*cos(d_psi)) + d_w*(cos(d_phi)*sin(d_theta)*sin(d_psi) - sin(d_phi)*cos(d_psi));
% dydt(12) = -d_u*sin(d_theta) + d_v*sin(d_phi)*cos(d_theta) + d_w*cos(d_phi)*cos(d_theta);


dydt = dydt';
