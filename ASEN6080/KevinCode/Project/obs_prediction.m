function yhat = obs_prediction(x_sc,r_obs)
    omegaE = 7.2921158553e-5;
    % calculate dRi in ECI
    dReci = x_sc(1:3) - r_obs(1:3)';
    
    % station velocties
    dxs = -r_obs(2)*omegaE;
    dys = r_obs(1)*omegaE;
    dzs = 0;
    Vs = [dxs dys dzs];
    
    % calculate dVi in ECI
    dVeci = x_sc(4:6) - Vs(1:3);

    % calculate rhoi from ECI frame
    rhoeci = norm(dReci);

    % range rate
    drhoeci = (dReci(1)*dVeci(1)+dReci(2)*dVeci(2)+dReci(3)*dVeci(3))/rhoeci;

    % predicted measurement
    yhat = [rhoeci; drhoeci];
end