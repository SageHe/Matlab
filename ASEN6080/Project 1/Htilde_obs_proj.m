function Hmat_obs = Htilde_obs_proj(x_sc,x_obs)
We = 7.2921158553e-5; % rotation rate of Earth in rads/s


x = x_sc(1);
y = x_sc(2);
z = x_sc(3);
xdot = x_sc(4);
ydot = x_sc(5);
zdot = x_sc(6);

xs = x_obs(1);
ys = x_obs(2);
zs = x_obs(3);
xsdot = x_obs(4);
ysdot = x_obs(5);
zsdot = x_obs(6);

Hmat_obs(1,:) = [-(abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), -(abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), -(abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];
Hmat_obs(2,:) = [-(abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), -(abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), -(abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];
end