function Htilde_sc_obs = Htilde_sc_rho_rhod(x_sc,x_obs)

    % assign variables
    % spacecraft state
    x = x_sc(1);
    y = x_sc(2);
    z = x_sc(3);
    xdot = x_sc(4);
    ydot = x_sc(5);
    zdot = x_sc(6);
    
    % observer station state
    xs = x_obs(1);
    ys = x_obs(2);
    zs = x_obs(3);
    dxs = x_obs(4);
    dys = x_obs(5);
    dzs = x_obs(6);
    
    % measurement sensitivity matrix
    Htilde_sc_obs = [                                                                                                                                             (abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                                              (abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                                              (abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                              0,                                                              0,                                                              0;
                    (abs(x - xs)*sign(x - xs)*((dxs - xdot)*(x - xs) + (dys - ydot)*(y - ys) + (dzs - zdot)*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (dxs - xdot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (abs(y - ys)*sign(y - ys)*((dxs - xdot)*(x - xs) + (dys - ydot)*(y - ys) + (dzs - zdot)*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (dys - ydot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (abs(z - zs)*sign(z - zs)*((dxs - xdot)*(x - xs) + (dys - ydot)*(y - ys) + (dzs - zdot)*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (dzs - zdot)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];

end