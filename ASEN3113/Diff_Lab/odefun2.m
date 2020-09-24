function yout = odefun2(t,y)

yout = zeros(2,1);

yout(1) = y(2);

mu = 8.53;
omega = ((2*pi)/10);
A = 1.2;

yout(2) = mu*(1 - y(1)^2)*y(2) - y(1) + A*sin(omega*t);
end