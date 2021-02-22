function [LLA] = ECEF_to_LLA(ECEF)

    lambda = atan2d(ECEF(2),ECEF(1));       % longitude, degrees
    rho = sqrt(ECEF(1)^2 + ECEF(2)^2);      % meters
    r = norm([ECEF(1) ECEF(2) ECEF(3)]);    % radius of ECEF coordinates, meters
    phiGD = asind(ECEF(3)/r);               % initial guess, degrees

    R = 6378.137*1000;      % meters
    f = 1/298.257223563;    % flattening parameter of ellipsoid
    e = sqrt(2*f-f^2);      % eccentricity of ellipsoid
    tol = 1;
    i = 1;
    
    while tol > 1e-10
        C = R/(1 - e^2*sind(phiGD(i))^2)^.5;                      % radius of curvature in the meridian, meters
        phiGD(i+1) = atan2d(ECEF(3) + C*e^2*sind(phiGD(i)),rho);  % update phiGD, degrees
        tol = phiGD(i+1) - phiGD(i);                              % update tolerance
        i = i+1;
    end
    
    phiGD = phiGD(i);
    h = rho/cosd(phiGD) - C; % ellipsoidal height, meters
    LLA = [phiGD lambda h];
end