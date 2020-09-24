function  dydt = ThermoModel(t,k)

%import global variables
global Cd pwater AThroat P0 Volair_i gamma P_amb ls g Theta p_air_amb CD ABottle Volbottle m_air_i R 

% set variables in diff eq vector
m_r = k(1); %Mass of rocket
m_air = k(2); %Mass of air
Vx = k(3); %Velocity in x direction
Vz = k(4); %Velocity in z direction
Vy = k(5); %Velocity in y direction
Px = k(6); %Position in x direction
Pz = k(7); %Position in z direction
Py = k(8); %Position in y direction
v_air = k(9); %Volume of air

% stop it once the rocket returns to the ground
if Pz<0
    dydt=[0 0 0 0 0 0 0 0 0];
    return
end


%Drag is universal; define it outside of if loops
V = sqrt((Vx^2)+(Vz^2)+(Vy^2));
D = (p_air_amb/2)*(V^2)*CD*ABottle;

%define whether it is on the rail or not with heading
if Px <= ls*cos(Theta)
    hx = cos(Theta);
    hz = sin(Theta);
    hy = 0;
else
    hx = Vx/V;
    hz = Vz/V;
    hy = Vy/V;
end



%% PHASE ONE: Water Propulsion
if v_air < Volbottle %If the volume of air is less than total volume 
    P = P0*(Volair_i/v_air)^gamma;     
    Thrust = 2*Cd*(P-P_amb)*AThroat;
    %acceleration
    ax = (Thrust - D)*hx/m_r;
    az = (Thrust - D)*hz/m_r-g;
    ay = (Thrust - D)*hy/m_r;
    
    dydt(1) = -Cd*AThroat*sqrt((2*pwater*(P-P_amb))); 
    dydt(2) = 0; 
    dydt(3) = ax;
    dydt(4) = az; 
    dydt(5) = ay;
    dydt(6) = Vx; 
    dydt(7) = Vz; 
    dydt(8) = Vy; 
    dydt(9) = Cd*AThroat*sqrt((2*(P-P_amb))/pwater);
    
%% PHASE TWO/THREE: Pressurized Air Exhaustion & Ballistic Phase
else 
    dydt(9) = 0;
    
    %Define Ending Variables of Pressure and Temperature
    Pend = P0*(Volair_i/Volbottle)^gamma; 

    %Define Pressure, Density, and Temperature for phase two
    P = ((m_air/m_air_i)^gamma)*Pend; 
    p = (m_air/Volbottle); 
    T = P/(p*R); 
    %PHASE 2
    if P > P_amb
        Pcr = P*(2/(gamma+1))^(gamma/(gamma-1));
        %SUBPHASE 2a | Pcr is Greater Than Ambient Pressure | Choked Flow
        if Pcr > P_amb
            Pexit = Pcr; 
            Texit = (2/(gamma+1))*T; 
            Vexit = sqrt(gamma*R*Texit); 
            pexit = Pexit/(R*Texit); 
            
            
        %SUBPHASE 2b | Pcr is Less Than or Equal to Ambient Pressure | Free Flow
        elseif Pcr <= P_amb
            Pexit = P_amb;  
            Mexit = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1)); 
            Texit = T/(1+((gamma-1)/2)*Mexit^2); 
            pexit = P_amb/(R*Texit); 
            Vexit = Mexit*sqrt(gamma*R*Texit); 
            
        end

        %force calculations
        Thrust = m_air*Vexit+(Pexit-P_amb)*AThroat;      
        %computing acceleration
        ax = (Thrust - D)*hx/(m_r);
        az = (Thrust - D)*hz/(m_r)-g;
        ay = (Thrust - D)*hy/(m_r);
        %mass flow is equivalent for rocket and air
        dydt(1) = -Cd*pexit*AThroat*Vexit;
        dydt(2) = -Cd*pexit*AThroat*Vexit;
    
%PHASE 3 : Ballistics
    elseif P <= P_amb
        %define acceleration without thrust
        ax = (-D)*hx/(m_r);
        az = (-D)*hz/(m_r)-g;
        ay = (-D)*hy/(m_r);
        dydt(1) = 0;
        dydt(2) = 0;
            
    end
    %these four are the same for both 2 and 3 phase
    dydt(3) = ax; 
    dydt(4) = az; 
    dydt(5) = ay;
    dydt(6) = Vx;
    dydt(7) = Vz;
    dydt(8) = Vy;
    
end
%make conditions for if it is in air or on ground 
if Pz > 0
    dydt(6) = Vx;
    dydt(7) = Vz;
    dydt(8) = Vy;
else
    dydt(6) = 0;
    dydt(7) = 0;
    dydt(8) = 0;
end

%transpose dydt matrix
dydt = dydt';

end


