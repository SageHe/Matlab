function zd = keplermuJ2J3_wPhi_ODE(t,z)

    mu = 3.986004415E5; % km^3/s^2
    
    % assign variables
    x = z(1);
    y = z(2);
    zz = z(3);
    xdot = z(4);
    ydot = z(5);
    zdot = z(6);    
    
    Re  = 6378.0; % km
    J2 = 0.0010826269;
    J3 = -0.0000025323;
    n = 9;
    
    % radius of s/c from center of Earth
    r = (x^2 + y^2 + zz^2)^.5;

    xddot_u = -mu*x/r^3;
    yddot_u = -mu*y/r^3;
    zddot_u = -mu*zz/r^3;

    xddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(zz/r)^2)*x/r);
    yddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(zz/r)^2)*y/r);
    zddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((3 - 5*(zz/r)^2)*zz/r);

    xddot_J3 = 1/2*J3*(mu/r^2)*(Re/r)^3*(5*(7*(zz/r)^3 - 3*(zz/r))*x/r);
    yddot_J3 = 1/2*J3*(mu/r^2)*(Re/r)^3*(5*(7*(zz/r)^3 - 3*(zz/r))*y/r);
    zddot_J3 = 1/2*J3*(mu/r^2)*(Re/r)^3*(3*(1 - 10*(zz/r)^2 + 35/3*(zz/r)^4));

    xddot = xddot_u + xddot_J2 + xddot_J3;
    yddot = yddot_u + yddot_J2 + yddot_J3;
    zddot = zddot_u + zddot_J2 + zddot_J3;
    
    % A matrix
    A = dfdZ_wJ2J3(z);
    
    % reshape STM
    Phi = reshape(z(n+1:end),n,n);

    % DKE
    Phi_dot = A*Phi;

    % reshape into (n^+n)x1 array to be integrated
    zd = [xdot; ydot; zdot; xddot; yddot; zddot; 0; 0; 0; reshape(Phi_dot,n^2,1)];
end