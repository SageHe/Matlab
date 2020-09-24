function [dydt] = scode(t,y)
theta = y(1);
omega = y(2);
theta_d = deg2rad(180);
domegadt = ((-4.775e-4)*(theta_d - theta) - .00219*omega)/(.0025);
dthetadt = 0;

dydt(1) = dthetadt;
dydt(2) = domegadt;

dydt = dydt';
end