function [Ab,Vb] = burn_geometry(r,h,rb)
    
    if rb >= r % motor is burnt out
        Ab = 0; % [m^2] 
    else % there is grain remaining
        %% BURN AREA
        Ab = ; % [m^2] total burn area

        %% BURN VOLUME
        % Vb = ; % [m^3]
    end