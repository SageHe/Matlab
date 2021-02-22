function [AZ, EL, RANGE] = azelrange(userECEF,userLLA,satECEF)

    LOS_ECEF = satECEF - userECEF;
    RANGE = norm(LOS_ECEF);
   
    %LOS_LLA = ECEF_to_LLA(userECEF);   % this function assumes an ellipsoidal Earth
    C_LOS_ENU_ECEF = ECEF_to_ENU(userLLA);
    LOS_ENU = C_LOS_ENU_ECEF*LOS_ECEF';
    
    AZ = atan2d(LOS_ENU(1),LOS_ENU(2));
    EL = asind(LOS_ENU(3)/norm(LOS_ENU));
    
end