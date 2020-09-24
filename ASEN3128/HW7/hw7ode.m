function [dydt] = hw7ode(t,y)
global A theta_0 u_0
delta_u = y(1);
delta_w = y(2);
delta_q = y(3);
delta_theta = y(4);
delta_x = y(5);
delta_z = y(6);

deltau_dot = A(1,1)*delta_u + A(1,2)*delta_w + A(1,4)*delta_theta;
deltaw_dot = A(2,1)*delta_u + A(2,2)*delta_w + A(2,3)*delta_q + A(2,4)*delta_theta;
deltaq_dot = A(3,1)*delta_u + A(3,2)*delta_w + A(3,3)*delta_q + A(3,4)*delta_theta;
deltatheta_dot = A(4,1)*delta_u + A(4,2)*delta_w + A(4,3)*delta_q + A(4,4)*delta_theta;
deltax_dot = delta_u*cos(theta_0) + delta_w*sin(theta_0) - u_0*delta_theta*sin(theta_0);
deltaz_dot = -delta_u*sin(theta_0) + delta_w*cos(theta_0) - u_0*delta_theta*cos(theta_0);

dydt(1) = deltau_dot;
dydt(2) = deltaw_dot;
dydt(3) = deltaq_dot;
dydt(4) = deltatheta_dot;
dydt(5) = deltax_dot;
dydt(6) = deltaz_dot;

dydt = dydt';
end