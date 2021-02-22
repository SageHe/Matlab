function Htilde_obs = Htilde_obs_rho_rhod(x_sc,x_obs)
    omegaE = 7.2921158553e-5;
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
    
    % measurement sensitivity matrix
    Htilde_obs = [
                                                                                                                                                                            -(abs(x - xs)*sign(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                                                                         -(abs(y - ys)*sign(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2),                                                                                                                                         -(abs(z - zs)*sign(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2);
                (abs(x - xs)*sign(x - xs)*((xdot + omegaE*ys)*(x - xs) + (ydot - omegaE*xs)*(y - ys) + zdot*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (xdot + omegaE*ys + omegaE*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2), (omegaE*xs - ydot + omegaE*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2) + (abs(y - ys)*sign(y - ys)*((xdot + omegaE*ys)*(x - xs) + (ydot - omegaE*xs)*(y - ys) + zdot*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2), (abs(z - zs)*sign(z - zs)*((xdot + omegaE*ys)*(x - xs) + (ydot - omegaE*xs)*(y - ys) + zdot*(z - zs)))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - zdot/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(1/2)];
            
end