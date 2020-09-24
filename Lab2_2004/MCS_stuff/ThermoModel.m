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
Dx = (p_air_amb/2)*(Vx^2)*CD*ABottle;
Dz = (p_air_amb/2)*(Vz^2)*CD*ABottle;

%define whether it is on the rail or not with heading
if Px <= ls*cos(Theta)
    hx = cos(Theta);
    hz = sin(Theta);
    Vy = 0;
else
    hx = Vx/sqrt((Vx^2)+(Vz^2)+(Vy^2));
    hz = Vz/sqrt((Vx^2)+(Vz^2)+(Vy)^2);
    hy = Vy/sqrt((Vx^2)+(Vz^2)+(Vy)^2);
end

% An alternative way to define drag
% V=sqrt(Vx^2+Vy^2+Vz^2);
% % calculate the drag force by using equation 2
% D=(p_air_amb/2)*(V^2)*CD*ABottle;;
% 
% %Drag is universal; define it outside of if loops
% Dx = D*hx;
% Dz = D*hz;


%% PHASE ONE: Water Propulsion
if v_air < Volbottle %If the volume of air is less than total volume 
    P = P0*(Volair_i/v_air)^gamma;     
    Thrustx = 2*Cd*(P-P_amb)*AThroat*hx;
    Thrustz = 2*Cd*(P-P_amb)*AThroat*hz;%EQ 8 
    %acceleration
    ax = (Thrustx - Dx)/m_r;
    az = ((Thrustz - Dz)/m_r)-g;
    
    dydt(1) = -Cd*AThroat*sqrt((2*pwater*(P-P_amb))); %EQ 4
    dydt(2) = 0; %mass of air is not changing
    dydt(3) = ax; %EQ 1
    dydt(4) = az; %EQ 1
    dydt(6) = Vx; %Integrated Velocity
    dydt(7) = Vz; %Integrated Velocity
    dydt(8) = Vy; %Integrated Velocity
    dydt(9) = Cd*AThroat*sqrt((2*(P-P_amb))/pwater); %OTHER EQ 10
    
%% PHASE TWO/THREE: Pressurized Air Exhaustion & Ballistic Phase
else %volume is equal to the bottle
    dydt(9) = 0;
    
    %Define Ending Variables of Pressure and Temperature
    Pend = P0*(Volair_i/Volbottle)^gamma; %EQ 13 
    %Tend = Tair_i*(Volair_i/Volbottle)^(gamma-1); %EQ 13

    %Define Pressure, Density, and Temperature for phase two
    P = ((m_air/m_air_i)^gamma)*Pend; %EQ 14
    p = (m_air/Volbottle); %EQ 15
    T = P/(p*R); %EQ 15
    %PHASE 2
    if P > P_amb
        Pcr = P*(2/(gamma+1))^(gamma/(gamma-1));
        %SUBPHASE 2a | Pcr is Greater Than Ambient Pressure | Choked Flow
        if Pcr > P_amb
            Pexit = Pcr; %Phase identity
            Texit = (2/(gamma+1))*T; %EQ 18
            Vexit = sqrt(gamma*R*Texit); %EQ 17
            pexit = Pexit/(R*Texit); %EQ 18
            
            
        %SUBPHASE 2b | Pcr is Less Than or Equal to Ambient Pressure | Free Flow
        elseif Pcr <= P_amb
            Pexit = P_amb; %Phase identity
            %Mexit = sqrt(2*((P/P_amb)^((gamma-1)/gamma)-1)/(gamma-1)); %EQ 20
            Mexit = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1)); %EQ 20 Andy's
            Texit = T/(1+((gamma-1)/2)*Mexit^2); %EQ 21
            pexit = P_amb/(R*Texit); %EQ 21
            Vexit = Mexit*sqrt(gamma*R*Texit); %EQ 22
            
        end
        %Exit velocity components
        Vexitx = Vexit*hx;
        Vexitz = Vexit*hz;
        %Fx = Cd*pexit*AThroat*Vexitx^2+(Pexit-P_amb)*AThroat;
        Thrustx = m_air*Vexitx+(Pexit-P_amb)*AThroat;
        Thrustz = m_air*Vexitz+(Pexit-P_amb)*AThroat;
        %Fz = Cd*pexit*AThroat*Vexitz^2+(Pexit-P_amb)*AThroat;        
        %computing acceleration
        ax = (Thrustx - Dx)/(m_r);
        az = (Thrustz - Dz)/(m_r)-g;
        %mass flow is equivalent for rocket and air
        dydt(1) = -Cd*pexit*AThroat*Vexit;
        dydt(2) = -Cd*pexit*AThroat*Vexit;
    
%PHASE 3 : Ballistics
    elseif P <= P_amb
        %define acceleration without thrust
        ax = (-Dx)/(m_r);
        az = (-Dz)/(m_r)-g;
        dydt(1) = 0;
        dydt(2) = 0;
            
    end
    %these four are the same for both 2 and 3 phase
    dydt(3) = ax; %EQ 1
    dydt(4) = az; %EQ 1
    dydt(6) = Vx;
    dydt(7) = Vz;
    dydt(8) = Vy;
    
end
%make conditions for if it is in air or on ground so that graph looks
%prettier
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


