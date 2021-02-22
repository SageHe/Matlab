function keps = car2kep(X,mu)
    % INPUTS
    % X:    6x1 or 1x6 state vector, units of SMA depend on input units
    % mu:   gravitational constant of central body, units consistent with state vector
    
    % OUTPUTS
    % keps: 1x6 vector [SMA ECC INC RAAN AP TA] (all angles output in degrees)

    r = [X(1); X(2); X(3)];
    v = [X(4); X(5); X(6)];
    
    h = cross(r,v);
    hmag = norm(h);
    rmag = norm(r);
    vmag = norm(v);
    ecc = cross(v,h)/mu - r/rmag;
    emag = norm(ecc);
    eps = .5*vmag^2 - mu/rmag;
    SMA = -mu/(2*eps);
    inc = acosd(h(3)/hmag);
    n = cross([0,0,1],h);
    nmag = norm(n);
    if dot(n,[0,1,0]) > 0
        RAAN = acosd(dot(n,[1,0,0])/nmag);
    else
        RAAN = -acosd(dot(n,[1,0,0])/nmag);
    end
    if dot(ecc,[0,0,1]) > 0
        ap = acosd(dot(n,ecc)/(nmag*emag));
    else
        ap = -acosd(dot(n,ecc)/(nmag*emag));
    end
    if dot(r,v) > 0
        TA = acosd(dot(r,ecc)/(rmag*emag));
    else
        TA = -acosd(dot(r,ecc)/(rmag*emag));
    end
    
    keps = [SMA emag inc RAAN ap TA];
    

end