
%Housekeeping
clc; clear; close all;
x_el = [];
y_el = [];
%% Set up given variables from input document
%Make all variables global
global mw mi mf Py0 windy windx windz g Cd p_air_amb Volbottle P_amb gamma pwater DThroat DBottle R MBottle CD P_gage Volwater_i Tair_i V0 Theta x0 z0 ls AThroat P0 Volair_i m_air_i ABottle counter
drag_coeff = linspace(.4,.6,10);
for i = 1:length(drag_coeff)
%All given variables, can be modified
g = 9.81;           %m/s2 ... acceleration due to gravity
Cd = 0.8;           %... discharge coefficient
p_air_amb = 0.961;  %kg/m^3 ... ambient air density
Volbottle = 0.002;  %m^3 ... volume of empty bottle
P_amb=12.1*6894.75; %pa ... atmospheric pressure
gamma = 1.4;        %... ratio of specific heats for air
pwater = 1000;      %kg/m^3 ... density of water
DThroat = 2.1/100;  %mm ... diameter of throat
DBottle = 10.5/100; %mm ... diameter of bottle
R = 287;            %J/kgK ... gas constant of air
P_gage = 40*6894.75;%pa ... initial gage pressure of air in bottle
V0 = 0.0;           %m/s ... initial velocity of rocket
Theta = 40*(pi/180);%radians ... initial angle of rocket
x0 = 0.0;           %m ... initial horizontal distance
z0 = 0.01;          %m ... initial vertical height
Py0 = 0;             %m ... initial sideways offset
ls= 0.5;            %m ... length of test stand

% The stuff that changes 
Tair_i = 278.9;       %K ... initial temperature of air
CD = drag_coeff(i);           %... drag coefficient
Volwater_i = 0.001; %m^3 ... initial volume of water inside bottle
mf = 0.124;         %final burnout mass of the rocket
MBottle = 0.124;     %kg ... mass of empty 2-liter bottle with cone and fins
mw = Volwater_i*pwater; % initial mass of the water


% The wind velocities in each direction
windx = 1.78816;
windy = 1;
windz = 0;

%% Calculate Derived Variables
%Initial Mass of Air
m_air_i = ((P_gage+P_amb)*(Volbottle-Volwater_i))/((R)*(Tair_i)); 
% initial mass of the rocket
mi = m_air_i + mw + mf;    
%Initial Total Pressure
P0 = P_gage+P_amb;
%Area of the Throat
AThroat = ((DThroat^2)*pi)/4; %m^2 ... 
%Volume of Air
Volair_i = Volbottle-Volwater_i; %m^3 ... Initial Volume of Air
%Cross-Sectional Area of the Bottle
ABottle = (pi*(DBottle^2))/4; %m^2 ... Area of Bottle

% Account for crosswinds and initial instantaneous deltaV
Vrelx = V0*cos(Theta)+windx;
Vrelz = V0*sin(Theta)+windz;
Vrely = windy;

%% Initial Values
%Initial value of total mass
y0(1) = mi;
%Initial value of air mass
y0(2) = m_air_i;
%Initial value of x-direction velocity
y0(3) = Vrelx;
%Initial value of z-direction velocity
y0(4) = Vrelz;
%Initial value of y-direction velocity
y0(5) = Vrely;
%Initial value of x-direction position
y0(6) = x0;
%Initial value of z-direction position
y0(7) = z0;
%Initial value of y-direction position
y0(8) = Py0;
%Initial value of volume of air with time
y0(9) = Volair_i;


%time span
tspan = [0 5];

%Call ODE45
[t,y] = ode45('ThermoModel',tspan, y0);

%Split up the results of ODE45
mr = y(:,1);
m_air = y(:,2);
Vx = y(:,3);
Vz = y(:,4);
Vy = y(:,5);
x = y(:,6);
z = y(:,7);
ypos = y(:,8)*2;
v = y(:,9);

% %Develop an array for thrust using two indexes
% Thrust = zeros(length(v),1);
% index1 = 0;
% index2 = 0;
% for i = 1:length(v)
%     if v(i) < Volbottle
%         P = P0*((Volair_i/v(i))^gamma);
%         Vexit = sqrt((2*(P-P_amb))/pwater);
%         massflow = Cd*pwater*AThroat*Vexit;
%         Thrust(i) = massflow*Vexit;
%         index1 = index1+1;
%         index2 = index2+1;
%     else
%         Pend = P0*(Volair_i/Volbottle)^gamma;
%         P = Pend * (m_air(i)/ m_air_i)^gamma;
%         p = m_air(i)/Volbottle;
%         T = P / (p*R);
%         if P > P_amb
%             Pcr = P* (2/(gamma+1))^(gamma/(gamma-1));
%             if Pcr > P_amb
%                 Pexit = Pcr; %Phase identity
%                 Texit = (2/(gamma+1))*T; %EQ 18
%                 Vexit = sqrt(gamma*R*Texit); %EQ 17
%                 pexit = Pcr/(R*Texit); %EQ 18
%             elseif Pcr <= P_amb
%                 Pexit = P_amb; %Phase identity
%                 %Mexit = sqrt(2*((P/P_amb)^((gamma-1)/gamma)-1)/(gamma-1)); %EQ 20
%                 Mexit = sqrt((2/(gamma-1))*((P/P_amb)^((gamma-1)/gamma)-1)); %EQ 20 Andy's
%                 Texit = T/(1+((gamma-1)/2)*Mexit^2); %EQ 21
%                 pexit = P_amb/(R*Texit); %EQ 21
%                 Vexit = Mexit*sqrt(gamma*R*Texit); %EQ 22
%             end
%             Thrust(i) = Vexit* m_air(i) + (Pexit-P_amb)*AThroat;
%             index2 = index2 + 1;
%         else
%             Thrust(i) = 0;
%         end
%     end
% end
% 
% 
% %Plot everything and make it look nice
% figure
% plot(t(1:index1),Thrust(1:index1),'c-')
% title('Thrust Profile')
% xlabel('Time (sec)')
% ylabel('Thrust (N)')
% hold on;
% plot(t(index1:index2),Thrust(index1:index2),'m-')
% hold on;
% plot(t(index2:(index2+10)),Thrust(index2:(index2+10)),'g-')
% legend('Phase One: Water Exhaustion','Phase Two: Air Pressure Exhaustion','Phase Three: Ballistic Flight')
% axis([0 .3 0 200])           
% 

% % find the Isp
% boindex = find(mr < 0.1260);
% boi = boindex(1);
% boVx = Vx(boi);
% boVz = Vz(boi);
% boV = sqrt(boVx^2 + boVz^2);
% Ispthermo = boV / (9.81*log(mi/mf));


%Again, plot everything
% figure
plot3(x,ypos,z)
hold on
axis([0 70 0 20])
title('Thermo Model for TA Baseline')
xlabel('X Distance (m)')
ylabel('Y Distance (m)')
zlabel('Altitude (m)')
fprintf('Bottle rocket''s maximum height is %f meters.\n\n Bottle rocket''s maximum distance is %f meters.\n',max(z),max(x));
grid on
x_el = [x_el x(end)];
y_el = [y_el ypos(end)];
end
% The TA flight went 48.829 meters. 
% The drift angle was 13 degrees.
realx = cos(deg2rad(13))*48.829;
realy = sin(deg2rad(13))*48.829;

modelx = max(x);
modely = max(ypos);

x_shift = mean(x_el);
y_shift = mean(y_el);
std_x = std(x_el);
std_y = std(y_el);
a = std_x; % horizontal radius
b = std_y; % vertical radius
x_center = x_shift; % x0,y0 ellipse centre coordinates
y_center = y_shift;
t=-pi:0.01:pi;
x=x_center + a*cos(t);
y=y_center + b*sin(t);
plot(x,y)

std_x = std(x_el);
std_y = std(y_el);
a = 2*std_x; % horizontal radius
b = 2*std_y; % vertical radius
x_center = x_shift; % x0,y0 ellipse centre coordinates
y_center = y_shift;
t = -pi:0.01:pi;
x = x_center + a*cos(t);
y = y_center + b*sin(t);
plot(x,y)
zlim([0 25])