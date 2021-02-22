function C_ENU_ECEF = ECEF_to_ENU(LLA)   

    % transformation matrix to convert ECEF to ENU
    C_ENU_ECEF = [-sind(LLA(2))                   cosd(LLA(2))                   0;
            -sind(LLA(1))*cosd(LLA(2))    -sind(LLA(1))*sind(LLA(2))    cosd(LLA(1));
            cosd(LLA(1))*cosd(LLA(2))     cosd(LLA(1))*sind(LLA(2))     sind(LLA(1))];
        
end