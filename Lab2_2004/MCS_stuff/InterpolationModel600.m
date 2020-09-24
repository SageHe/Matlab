function  dk = InterpolationModel600(t,j)

%import global variables
global mthermo mtimes F mf g CD p_air_amb Theta z0 ls ABottle

%% set up the equations to be solved

x = j(1); % x distance
y = j(2); % y distance
z = j(3); % z distance
Vx = j(4); % the velocity of x axis
Vy = j(5); % the velocity of y axis
Vz = j(6); % the velocity of z axis

%% once the height is smaller than 0, the water bottle rocket stops moving
% and all the change of variables equal to zero
if z<0
    dk=[0 0 0 0 0 0];
    return
end

%% calculate the total Velocity
V=sqrt(Vx^2+Vy^2+Vz^2);

%% Calculate Drag
% calculate the drag force by using equation 2
D=0.5*p_air_amb*V^2*CD*ABottle;

%% Check if the rocket is on the test stand
if sqrt(x^2+(z-z0)^2)<ls || sqrt(x^2+(z-z0)^2)==ls % test stand length
    % if the bottle rocket is on the test stand   
    hx=cos(Theta);
    hz=sin(Theta);
    hy=0;
else
    % if the bottle rocket leaves the test stand
    hx=Vx/V; % cos theta
    hz=Vz/V; % sin theta
    hy=Vy/V;
end


%% Determine the Forces Acting on the Rocket
%if t<tThrust
    
% call thrust interpolation
T = ThrustInterpolation600(t,F);

% Find the Mass of the Rocket
m = MassInterpolation(t,mthermo,mtimes,mf);

% Define the accelerations over the flight
ax = (T - D)*hx/m;
ay = (T - D)*hy/m;
az = (T - D)*hz/m - g;
    
% else
%     % thrust cuts out
%     T = 0;
%     % mass of the rocket does not change anymore
%     m = mf;
%     % Define the accelerations over the flight
%     ax = (-D)*hx/m;
%     ay = (-D)*hy/m;
%     az = (-D)*hz/m - g;
% end


%% Write the equations to be solved

dk(1) = Vx;  
dk(2) = Vy;  
dk(3) = Vz; 
dk(4) = ax; 
dk(5) = ay; 
dk(6) = az; 

dk=dk';

%[Vx,Vy,Vz,m,t]

end



   