function dfdz = dfdZ_wJ2J3(Z)
    % assign variables
    x = Z(1);
    y = Z(2);
    z = Z(3);
    mu = Z(7);
    J2 = Z(8);
    J3 = Z(9);
    
    Re  = 6378.0; % km
    
    % intermediate
    r = (x^2 + y^2 + z^2)^(1/2);

    % partials
    dfdz = [                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 1, 0, 0,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 1, 0,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 1,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0;
            (3*mu*x^2)/r^5 - mu/r^3 - (J3*Re^3*mu*((15*z)/r - (35*z^3)/r^3))/(2*r^6) + (3*J2*Re^2*mu*((5*z^2)/r^2 - 1))/(2*r^5) + (J3*Re^3*mu*x*((15*x*z)/r^3 - (105*x*z^3)/r^5))/(2*r^6) - (15*J2*Re^2*mu*x^2*z^2)/r^9 - (15*J2*Re^2*mu*x^2*((5*z^2)/r^2 - 1))/(2*r^7) + (3*J3*Re^3*mu*x^2*((15*z)/r - (35*z^3)/r^3))/r^8,                                                                                                                                                                                                                      (3*mu*x*y)/r^5 + (J3*Re^3*mu*x*((15*y*z)/r^3 - (105*y*z^3)/r^5))/(2*r^6) + (3*J3*Re^3*mu*x*y*((15*z)/r - (35*z^3)/r^3))/r^8 - (15*J2*Re^2*mu*x*y*z^2)/r^9 - (15*J2*Re^2*mu*x*y*((5*z^2)/r^2 - 1))/(2*r^7),                                                                                                     (3*mu*x*z)/r^5 + (3*J2*Re^2*mu*x*((10*z)/r^2 - (10*z^3)/r^4))/(2*r^5) - (J3*Re^3*mu*x*(15/r - (120*z^2)/r^3 + (105*z^4)/r^5))/(2*r^6) + (3*J3*Re^3*mu*x*z*((15*z)/r - (35*z^3)/r^3))/r^8 - (15*J2*Re^2*mu*x*z*((5*z^2)/r^2 - 1))/(2*r^7), 0, 0, 0, (3*J2*Re^2*x*((5*z^2)/r^2 - 1))/(2*r^5) - (J3*Re^3*x*((15*z)/r - (35*z^3)/r^3))/(2*r^6) - x/r^3, (3*Re^2*mu*x*((5*z^2)/r^2 - 1))/(2*r^5), -(Re^3*mu*x*((15*z)/r - (35*z^3)/r^3))/(2*r^6);
                                                                                                                                                                                                                     (3*mu*x*y)/r^5 + (J3*Re^3*mu*y*((15*x*z)/r^3 - (105*x*z^3)/r^5))/(2*r^6) + (3*J3*Re^3*mu*x*y*((15*z)/r - (35*z^3)/r^3))/r^8 - (15*J2*Re^2*mu*x*y*z^2)/r^9 - (15*J2*Re^2*mu*x*y*((5*z^2)/r^2 - 1))/(2*r^7), (3*mu*y^2)/r^5 - mu/r^3 - (J3*Re^3*mu*((15*z)/r - (35*z^3)/r^3))/(2*r^6) + (3*J2*Re^2*mu*((5*z^2)/r^2 - 1))/(2*r^5) + (J3*Re^3*mu*y*((15*y*z)/r^3 - (105*y*z^3)/r^5))/(2*r^6) - (15*J2*Re^2*mu*y^2*z^2)/r^9 - (15*J2*Re^2*mu*y^2*((5*z^2)/r^2 - 1))/(2*r^7) + (3*J3*Re^3*mu*y^2*((15*z)/r - (35*z^3)/r^3))/r^8,                                                                                                     (3*mu*y*z)/r^5 + (3*J2*Re^2*mu*y*((10*z)/r^2 - (10*z^3)/r^4))/(2*r^5) - (J3*Re^3*mu*y*(15/r - (120*z^2)/r^3 + (105*z^4)/r^5))/(2*r^6) + (3*J3*Re^3*mu*y*z*((15*z)/r - (35*z^3)/r^3))/r^8 - (15*J2*Re^2*mu*y*z*((5*z^2)/r^2 - 1))/(2*r^7), 0, 0, 0, (3*J2*Re^2*y*((5*z^2)/r^2 - 1))/(2*r^5) - (J3*Re^3*y*((15*z)/r - (35*z^3)/r^3))/(2*r^6) - y/r^3, (3*Re^2*mu*y*((5*z^2)/r^2 - 1))/(2*r^5), -(Re^3*mu*y*((15*z)/r - (35*z^3)/r^3))/(2*r^6);
                                                                                                                                                                                                                         (3*mu*x*z)/r^5 + (J3*Re^3*mu*((60*x*z^2)/r^4 - (140*x*z^4)/r^6))/(2*r^5) - (15*J2*Re^2*mu*x*z^3)/r^9 - (5*J3*Re^3*mu*x*((35*z^4)/r^4 - (30*z^2)/r^2 + 3))/(2*r^7) - (15*J2*Re^2*mu*x*z*((5*z^2)/r^2 - 3))/(2*r^7),                                                                                                                                                                                                                          (3*mu*y*z)/r^5 + (J3*Re^3*mu*((60*y*z^2)/r^4 - (140*y*z^4)/r^6))/(2*r^5) - (15*J2*Re^2*mu*y*z^3)/r^9 - (5*J3*Re^3*mu*y*((35*z^4)/r^4 - (30*z^2)/r^2 + 3))/(2*r^7) - (15*J2*Re^2*mu*y*z*((5*z^2)/r^2 - 3))/(2*r^7), (3*mu*z^2)/r^5 - mu/r^3 - (J3*Re^3*mu*((60*z)/r^2 - (200*z^3)/r^4 + (140*z^5)/r^6))/(2*r^5) + (3*J2*Re^2*mu*((5*z^2)/r^2 - 3))/(2*r^5) + (3*J2*Re^2*mu*z*((10*z)/r^2 - (10*z^3)/r^4))/(2*r^5) - (15*J2*Re^2*mu*z^2*((5*z^2)/r^2 - 3))/(2*r^7) - (5*J3*Re^3*mu*z*((35*z^4)/r^4 - (30*z^2)/r^2 + 3))/(2*r^7), 0, 0, 0,   (J3*Re^3*((35*z^4)/r^4 - (30*z^2)/r^2 + 3))/(2*r^5) - z/r^3 + (3*J2*Re^2*z*((5*z^2)/r^2 - 3))/(2*r^5), (3*Re^2*mu*z*((5*z^2)/r^2 - 3))/(2*r^5),    (Re^3*mu*((35*z^4)/r^4 - (30*z^2)/r^2 + 3))/(2*r^5);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                                                                                                                                                                                               0,                                                                         0,                                                                                                        0];
end