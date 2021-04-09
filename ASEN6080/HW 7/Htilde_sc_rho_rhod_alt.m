function Htilde_sc_rho_rhod = Htilde_sc_rho_rhod_alt(x_sc,x_obs)
x = x_sc(1);
y = x_sc(2);
z = x_sc(3);
xDot = x_sc(4);
yDot = x_sc(5);
zDot = x_sc(6);

xS = x_obs(1);
yS = x_obs(2);
zS = x_obs(3);
xSDot = x_obs(4);
ySDot = x_obs(5);
zSDot = x_obs(6);

% Htilde_sc_rho_rhod(1,:) = [(abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), 0, 0, 0];
% Htilde_sc_rho_rhod(2,:) = [(xdot - xsdot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2) - (abs(x - xs)*sign(x - xs)*((x - xs)*(xdot - xsdot) + (y - ys)*(ydot - ysdot) + (z - zs)*(zdot - zsdot)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2), (ydot - ysdot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2) - (abs(y - ys)*sign(y - ys)*((x - xs)*(xdot - xsdot) + (y - ys)*(ydot - ysdot) + (z - zs)*(zdot - zsdot)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2), (zdot - zsdot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2) - (abs(z - zs)*sign(z - zs)*((x - xs)*(xdot - xsdot) + (y - ys)*(ydot - ysdot) + (z - zs)*(zdot - zsdot)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2), (x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];
dH1 = [(2*x - 2*xS)/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2)), (2*y - 2*yS)/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2)), (2*z - 2*zS)/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2)), 0, 0, 0];
dH2 = [(xDot - xSDot)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2) - ((2*x - 2*xS)*((x - xS)*(xDot - xSDot) + (y - yS)*(yDot - ySDot) + (z - zS)*(zDot - zSDot)))/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(3/2)), (yDot - ySDot)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2) - ((2*y - 2*yS)*((x - xS)*(xDot - xSDot) + (y - yS)*(yDot - ySDot) + (z - zS)*(zDot - zSDot)))/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(3/2)), (zDot - zSDot)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2) - ((2*z - 2*zS)*((x - xS)*(xDot - xSDot) + (y - yS)*(yDot - ySDot) + (z - zS)*(zDot - zSDot)))/(2*((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(3/2)), (x - xS)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2), (y - yS)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2), (z - zS)/((x - xS)^2 + (y - yS)^2 + (z - zS)^2)^(1/2)];

Htilde_sc_rho_rhod(1,:) = dH1;
Htilde_sc_rho_rhod(2,:) = dH2;

end