% Declan Murray
% 104700338
% ASEN 3128 - HW 7

%% Housekeeping
clc, clear all, close all;

%% introduce calculated x,z,m values from part 1
Xu = -8468.02; Xw = -8669.1; Xq = -76893; Xww = 481.5;
Zu = -53092; Zw = -124012.6; Zq = -644841.5; Zww = 4038.4;
Mu = -3851.1; Mw = -225533.8; Mq = -2183820.9; Mww = -36254;
u0 = 157.886;% m/s
g = 9.81;%m/s
theta0 = 0;
% subscript ww denotes \dot(w)

%% Introduce Inertia and mass needed to calculate A matrix
m = 2.83e6/9.81;% kg - total weight of the aircraft divided by grav const.
Ix = 2.46e7;% kg-m^2
Iy = 4.48e7;%kg-m^2
Iz = 6.73e7;%kg-m^2
Izx = 13.2e5;%kg-m^2

%% construct the A matrix

A = zeros(4);
% going column by column for debugging purposes
A(:,1) = [Xu/m; Zu/(m - Zww); 1/Iy * (Mu + Mww*Zu/(m - Zww)); 0];
A(:,2) = [Xw/m; Zw/(m - Zww); 1/Iy * (Mw + Mww*Zw/(m - Zww)); 0];
A(:,3) = [0; (Zq + m*u0)/(m - Zww); 1/Iy * (Mq + Mww*(Zq+m*u0)/(m-Zww)); 1];
A(:,4) = [-g*cos(theta0); -m*g*sin(theta0)/(m - Zww); -Mww*m*g*sin(theta0)/(Iy*(m-Zww)); 0];

%% Find Eigenvectors and Eigenvals

[eVecs, eVals] = eig(A);

%% calculate modal frequency and damping ratio

wS = eVals(1,1);
wP = eVals(3,3);

wnS = norm(eVals(1,1));
wnP = norm(eVals(3,3));

zetaS = -real(wS)/wnS;
zetaP = -real(wP)/wnP;

%period of phugoid mode
Tp = 2*pi/wnP;

%% calculate approximation eigenvalues
lam1 = Mq/(2*Iy) + 1/(2*Iy) *sqrt(Mq^2 + 4*Iy*u0*Mw);
lam2 = Mq/(2*Iy) - 1/(2*Iy) *sqrt(Mq^2 + 4*Iy*u0*Mw);

%lancaster approx period
Tl = pi*sqrt(2) * u0/g;

%% QUESTION 5

%place X,z,q vals into single matrix for easy use
X = [Xu Xw Xq Xww];
Z = [Zu Zw Zq Zww];
M = [Mu Mw Mq Mww];

constVals = [u0 theta0 m Iy X Z M];

%introduce initial conditions
du = 0;
dw = 0;
dq = 0;
dtheta = 0.1;
dxdot = 0;
dzdot = 0;

y0 = [du dw dq dtheta dxdot dzdot];

%time span
tspan = [0 100];

% call ODE funciton
[t, y] = ode45(@(t,y) hw7_ODE(t,y,constVals), tspan, y0);

%% Plot results
figure
% U vel
subplot(3,2,1)
plot(t,y(:,1));
title('U Component of Velocity over Time')
xlabel('Time (s)')
ylabel('U velocity (m/s)')
% W vel
subplot(3,2,2)
plot(t,y(:,2));
title('W Component of Velocity over Time')
xlabel('Time (s)')
ylabel('W velocity (m/s)')
%q
subplot(3,2,3)
plot(t,y(:,3));
title('Pitch, q, over Time')
xlabel('Time (s)')
ylabel('q (rad/s)')
% theta
subplot(3,2,4)
plot(t,y(:,4));
title('Theta over Time')
xlabel('Time (s)')
ylabel('Theta Value (rad)')

% xloc
subplot(3,2,5)
plot(t,y(:,5));
title('X Position over Time')
xlabel('Time (s)')
ylabel('X position (m)')

% z loc
subplot(3,2,6)
plot(t,y(:,6));
title('Z Position over Time')
xlabel('Time (s)')
ylabel('Z position (m)')













