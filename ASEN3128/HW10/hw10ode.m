function [dydt] = hw10ode(t,y)
global A theta_0 u0
delta_v = y(1);
delta_p = y(2);
delta_r = y(3);
delta_phi = y(4);
delta_psi = y(5);
delta_y = y(6);

deltav_dot = A(1,1)*delta_v + A(1,2)*delta_p + A(1,3)*delta_r + A(1,4)*delta_phi;
deltap_dot = A(2,1)*delta_v + A(2,2)*delta_p + A(2,3)*delta_r + A(2,4)*delta_phi;
deltar_dot = A(3,1)*delta_v + A(3,2)*delta_p + A(3,3)*delta_r + A(3,4)*delta_phi;
deltaphi_dot = A(4,1)*delta_v + A(4,2)*delta_p + A(4,3)*delta_r + A(4,4)*delta_phi;
deltapsi_dot = delta_r*sec(theta_0);
deltay_dot = u0*delta_psi*cos(theta_0) + delta_v;

dydt(1) = deltav_dot;
dydt(2) = deltap_dot;
dydt(3) = deltar_dot;
dydt(4) = deltaphi_dot;
dydt(5) = deltapsi_dot;
dydt(6) = deltay_dot;

dydt = dydt';
end