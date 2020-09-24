function [Ab,Vb] = burn_geometry(r,h,rb,h_grain)
    
    if rb >= r % motor is burnt out
        Ab = 0; % [m^2] 
    else % there is grain remaining
        %% BURN AREA
        Ab = 2*pi*rb*h + 2*(pi*(r^2 - rb^2)); % [m^2] total burn area

        %% BURN VOLUME
        Vb = pi*r^2*h_grain - pi*(r^2 - rb^2)*h; % [m^3]
    end