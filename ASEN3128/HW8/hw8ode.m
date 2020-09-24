function [dydt] = hw8ode(t,y)
global A theta_0 u_0 mat
delta_u = y(1);
delta_w = y(2);
delta_q = y(3);
delta_theta = y(4);
delta_x = y(5);
delta_z = y(6);

Y_dot = mat*[delta_u delta_w delta_q delta_theta]';
deltax_dot = delta_u*cos(theta_0) + delta_w*sin(theta_0) - u_0*delta_theta*sin(theta_0);
deltaz_dot = -delta_u*sin(theta_0) + delta_w*cos(theta_0) - u_0*delta_theta*cos(theta_0);

dydt(1) = Y_dot(1);
dydt(2) = Y_dot(2);
dydt(3) = Y_dot(3);
dydt(4) = Y_dot(4);
dydt(5) = deltax_dot;
dydt(6) = deltaz_dot;

dydt = dydt';
end