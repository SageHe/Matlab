function dydt = odefun(t,y0,mu)
dydt = zeros(4,1);

x = y0(1);
y = y0(2);
xd = y0(3);
yd = y0(4);


dvdx = x + ((mu - 1)*(x+mu))/(((x+mu)^2 + y^2)^(3/2)) - (mu*(x - 1 + mu))/(((x - 1 + mu)^2 + y^2)^(3/2));
dvdy = y + ((mu-1)*y)/(((x + mu)^2 + y^2)^(3/2)) - (mu*y)/(((x - 1 + mu)^2 + y^2)^(3/2));

dydt(1) = xd;
dydt(2) = yd;
dydt(3) = 2*yd + dvdx;
dydt(4) = -2*xd + dvdy;
end