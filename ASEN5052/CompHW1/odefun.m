function dydt = odefun(t,y,mu)
dydt = zeros(6,1);
r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
dydt(4) = (-mu*y(1))/r^3;
dydt(5) = (-mu*y(2))/r^3;
dydt(6) = (-mu*y(3))/r^3;
end