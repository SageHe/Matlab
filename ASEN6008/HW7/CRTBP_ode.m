function dydt = CRTBP_ode(Z)
dydt = zeros(6,1);

mu = 0.012150585609624; 

x = Z(1);
y = Z(2);
z = Z(3);
xdot = Z(4);
ydot = Z(5);
zdot = Z(6);

R1 = sqrt((x + mu)^2 + y^2 + z^2);
R2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

xdd = 2*ydot + x - (1 - mu)*((x+mu)/(R1^3)) - mu*((x - 1 + mu)/(R2^3));
ydd = -2*xdot + y - (1 - mu)*(y/(R1^3)) - mu*(y/(R2^3));
zdd = -(1 - mu)*(z/(R1^3)) - mu*(z/(R2^3));

dydt = [xdot ydot zdot xdd ydd zdd]';
end
