function zd = keplerJ2_wPhi_ODE(t,z)

    mu = 3.986004415E5; % km^3/s^2
    
    % assign variables
    x = z(1);
    y = z(2);
    zz = z(3);
    xdot = z(4);
    ydot = z(5);
    zdot = z(6);    
    
    Re  = 6378.0; % km
    J2 = 1.082626925638815e-3; 
    n = 7;
    
    % radius of s/c from center of Earth
    r = (x^2 + y^2 + zz^2)^(1/2);
    
    % accelerations due to mass of Earth
    xddot_u = -mu*x/r^3;
    yddot_u = -mu*y/r^3;
    zddot_u = -mu*zz/r^3;
    
    % accelerations due to J2
    xddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(zz/r)^2)*x/r);
    yddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(zz/r)^2)*y/r);
    zddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((3 - 5*(zz/r)^2)*zz/r);
    
    % total accelerations
    xddot = xddot_u + xddot_J2;
    yddot = yddot_u + yddot_J2;
    zddot = zddot_u + zddot_J2;    
    
    % A matrix
    A = [                                                                                                                                                                                                                                                                               0,                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                       0, 1, 0, 0,                                                                        0;
                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                       0, 0, 1, 0,                                                                         0;
                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                       0, 0, 0, 1,                                                                         0;
    (3*mu*x^2)/r^5 - mu/r^3 + (3*J2*Re^2*mu*((5*zz^2)/r^2 - 1))/(2*r^5) - (15*J2*Re^2*mu*x^2*zz^2)/r^9 - (15*J2*Re^2*mu*x^2*((5*zz^2)/r^2 - 1))/(2*r^7),                                                                                                           (3*mu*x*y)/r^5 - (15*J2*Re^2*mu*x*y*zz^2)/r^9 - (15*J2*Re^2*mu*x*y*((5*zz^2)/r^2 - 1))/(2*r^7),                                                                                                           (3*mu*x*zz)/r^5 + (3*J2*Re^2*mu*x*((10*zz)/r^2 - (10*zz^3)/r^4))/(2*r^5) - (15*J2*Re^2*mu*x*zz*((5*zz^2)/r^2 - 1))/(2*r^7), 0, 0, 0, (3*Re^2*mu*x*((5*zz^2)/r^2 - 1))/(2*r^5);
                                                                                                              (3*mu*x*y)/r^5 - (15*J2*Re^2*mu*x*y*zz^2)/r^9 - (15*J2*Re^2*mu*x*y*((5*zz^2)/r^2 - 1))/(2*r^7), (3*mu*y^2)/r^5 - mu/r^3 + (3*J2*Re^2*mu*((5*zz^2)/r^2 - 1))/(2*r^5) - (15*J2*Re^2*mu*y^2*zz^2)/r^9 - (15*J2*Re^2*mu*y^2*((5*zz^2)/r^2 - 1))/(2*r^7),                                                                                                           (3*mu*y*zz)/r^5 + (3*J2*Re^2*mu*y*((10*zz)/r^2 - (10*zz^3)/r^4))/(2*r^5) - (15*J2*Re^2*mu*y*zz*((5*zz^2)/r^2 - 1))/(2*r^7), 0, 0, 0, (3*Re^2*mu*y*((5*zz^2)/r^2 - 1))/(2*r^5);
                                                                                                                (3*mu*x*zz)/r^5 - (15*J2*Re^2*mu*x*zz^3)/r^9 - (15*J2*Re^2*mu*x*zz*((5*zz^2)/r^2 - 3))/(2*r^7),                                                                                                             (3*mu*y*zz)/r^5 - (15*J2*Re^2*mu*y*zz^3)/r^9 - (15*J2*Re^2*mu*y*zz*((5*zz^2)/r^2 - 3))/(2*r^7), (3*mu*zz^2)/r^5 - mu/r^3 + (3*J2*Re^2*mu*((5*zz^2)/r^2 - 3))/(2*r^5) + (3*J2*Re^2*mu*zz*((10*zz)/r^2 - (10*zz^3)/r^4))/(2*r^5) - (15*J2*Re^2*mu*zz^2*((5*zz^2)/r^2 - 3))/(2*r^7), 0, 0, 0, (3*Re^2*mu*zz*((5*zz^2)/r^2 - 3))/(2*r^5);
                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                       0, 0, 0, 0,                                                                         0];
    % reshape STM
    Phi = reshape(z(n+1:end),n,n);

    % DKE
    Phi_dot = A*Phi;

    % reshape into (n^+n)x1 array to be integrated
    zd = [xdot; ydot; zdot; xddot; yddot; zddot; 0; reshape(Phi_dot,n^2,1)];
end