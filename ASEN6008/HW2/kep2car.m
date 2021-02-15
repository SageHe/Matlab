function x_car = kep2car(OEs, mu, angle_type)
    % this function transforms Keplerian orbital elements into Cartesian coordinates
    
    % INPUTS
    % OEs:          6x1 or 1x6 vector of the six orbital elements (a, e, i RAAN, arg of per, true anom)
    % mu:           mass parameter of central body
    % angle_type:   can take in degrees or radians
    
    % OUTPUTS
    % x_car:        1x6 vector of Cartesian coordinates [x, y, z, xdot, ydot, zdot];
    %                                                   [km,km,km,km/s, km/s, km/s]
    
    SMA = OEs(1);
    ecc = OEs(2);
    inc = OEs(3);
    RAAN= OEs(4);
    arg = OEs(5);
    TA  = OEs(6);
    
    if angle_type == 'degrees'
        tha = arg + TA;
        radius = (SMA*(1-ecc^2))/(1+ecc*cosd(TA));
        vel = sqrt(2*mu/radius - mu/SMA);
        FPA = ecc*sind(TA)/(1+ecc*cosd(TA))*360/2/pi;
        vr = vel*sind(FPA);
        vtha = vel*cosd(FPA);
        rx = radius*(cosd(RAAN)*cosd(tha) - sind(RAAN)*cosd(inc)*sind(tha));
        ry = radius*(sind(RAAN)*cosd(tha) + cosd(RAAN)*cosd(inc)*sind(tha));
        rz = radius*(sind(inc)*sind(tha));
        r_car = [rx ry rz];
        vx = vr*(cosd(RAAN)*cosd(tha) - sind(RAAN)*cosd(inc)*sind(tha)) + vtha*(-cosd(RAAN)*sind(tha)-sind(RAAN)*cosd(inc)*cosd(tha));
        vy = vr*(sind(RAAN)*cosd(tha) + cosd(RAAN)*cosd(inc)*sind(tha)) + vtha*(-sind(RAAN)*sind(tha)+cosd(RAAN)*cosd(inc)*cosd(tha));
        vz = vr*(sind(inc)*sind(tha)) + vtha*sind(inc)*cosd(tha);
        v_car = [vx vy vz];
        x_car = [r_car v_car];
    elseif angle_type == 'radians'
        tha = arg + TA;
        radius = (SMA*(1-ecc^2))/(1+ecc*cos(TA));
        vel = sqrt(2*mu/radius - mu/SMA);
        FPA = ecc*sin(TA)/(1+ecc*cos(TA));
        vr = vel*sin(FPA);
        vtha = vel*cos(FPA);
        rx = radius*(cos(RAAN)*cos(tha) - sin(RAAN)*cos(inc)*sin(tha));
        ry = radius*(sin(RAAN)*cos(tha) + cos(RAAN)*cos(inc)*sin(tha));
        rz = radius*(sin(inc)*sin(tha));
        r_car = [rx ry rz];
        vx = vr*(cos(RAAN)*cos(tha) - sin(RAAN)*cos(inc)*sin(tha)) + vtha*(-cos(RAAN)*sin(tha)-sin(RAAN)*cos(inc)*cos(tha));
        vy = vr*(sin(RAAN)*cos(tha) + cos(RAAN)*cos(inc)*sin(tha)) + vtha*(-sin(RAAN)*sin(tha)+cos(RAAN)*cos(inc)*cos(tha));
        vz = vr*(sin(inc)*sin(tha)) + vtha*sin(inc)*cos(tha);
        v_car = [vx vy vz];
        x_car = [r_car v_car];
    else
        disp('Please enter the either "degrees" or "radians" for angle_type')
    end
end