function dydt = odefunP6(t,y0,mu)
dydt = zeros(20,1);

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

% state = [x;y;xd;yd];
% EOM = [dydt(1);dydt(2);dydt(3);dydt(4)];

Amat = zeros(4,4);
Amat(1:2,3:4) = eye(2);
Amat(3,1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1;
Amat(4,1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
Amat(3,2) = (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2);
Amat(4,2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
Amat(3,4) = 2;
Amat(4,3) = -2;

stm = reshape(y0(5:end),4,4);

stmdot = Amat*stm;
dydt(5:end) = reshape(stmdot,16,1);

end