function yhat = obs_prediction(x_sc,x_obs)
    % x_sc, x_obs must be a 1x6 vectors in ECI

    % calculate dRi in ECI
    dReci = x_sc(1:3) - x_obs(1:3);
    
    % calculate rhoi from ECI frame
    rhoeci = norm(dReci);
   
    % calculate dVi in ECI
    dVeci = x_sc(4:6) - x_obs(4:6);

    % range rate
    drhoeci = (dReci(1)*dVeci(1)+dReci(2)*dVeci(2)+dReci(3)*dVeci(3))/rhoeci;

    % predicted measurement
    yhat = [rhoeci; drhoeci];
end