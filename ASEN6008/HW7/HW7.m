clear all;close all;clc
%% Problem 1
%convert state and time into dimensionless quantities.
mu_S = 1.32712440018e11; %km
mu_EM = 4.035032351966808e5; %km

Pos = [150000000 6000 1450]'; %km

Vel = [0.00075 0.08 0.019]'; %km/s

Time = 450*86400; %Days

AU = 1.49597870691e8; %km

P = 2*pi*sqrt((AU^3)/(mu_S + mu_EM));
n = (2*pi)/P;

Time_NonDim = Time*n;

TimeUnit = Time/Time_NonDim;

Pos_NonDim = Pos./AU;
Vel_NonDim = Vel.*(TimeUnit/AU);
%% Problem 2
%Integrate states for their respective prop. times and plot each traj. in
%rotating frame
DistanceUnit = 384747.962856037;
mu = 0.012150585609624;
ICs = [1.2 0 0 0 -1.06110124 0 6.20628;...
        0.85 0 0.17546505 0 0.2628980369 0 2.5543991];
opts = odeset('RelTol',10e-12,'AbsTol',10e-12);

    
tspan = [0 ICs(1,end)];
Z = ICs(1,1:6);
    
[t,state1] = ode45(@(t,Z) CRTBP_ode(Z),tspan,Z,opts);
Z = ICs(2,1:6);
tspan = [0 ICs(2,end)];
[t,state2] = ode45(@(t,Z) CRTBP_ode(Z),tspan,Z,opts);

%scale back up
state1(:,1:3) = state1(:,1:3)*DistanceUnit;
state2(:,1:3) = state2(:,1:3)*DistanceUnit;

figure
plot3(state1(:,1),state1(:,2),state1(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
title('Trajectory 1, Rotating Frame')
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
legend('Trajectory','Earth','Moon')

figure
plot3(state2(:,1),state2(:,2),state2(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
title('Trajectory 2, Rotating Frame')
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
legend('Trajectory','Earth','Moon')
%% Problem 3
% Int. the ICs in table 2 and plot trajs. in both rotating and
% Earth-centered inertial frames. 

ICs = [-0.083 -0.03 0.01 3.53 -3.1 -0.1 26;...
       0.05 -0.05 0 4.0 2.6 0 15;...
       0.875 0 0.1914 0 0.23234 0 15;...
       -0.05 -0.02 0 4.09 -5.27 0 2.5];
   
tspan1 = [0 ICs(1,end)];
Z1 = ICs(1,1:6);
[t1,state1] = ode45(@(t,Z) CRTBP_ode(Z),tspan1,Z1,opts);
tspan2 = [0 ICs(2,end)];
Z2 = ICs(2,1:6);
[t2,state2] = ode45(@(t,Z) CRTBP_ode(Z),tspan2,Z2,opts);
tspan3 = [0 ICs(3,end)];
Z3 = ICs(3,1:6);
[t3,state3] = ode45(@(t,Z) CRTBP_ode(Z),tspan3,Z3,opts);
tspan4 = [0 ICs(4,end)];
Z4 = ICs(4,1:6);
[t4,state4] = ode45(@(t,Z) CRTBP_ode(Z),tspan4,Z4,opts);

theta0 = 0;
theta_dot = 1;

state1(:,1) = state1(:,1) + mu;
state1(:,1:3) = state1(:,1:3)*DistanceUnit;
state2(:,1) = state2(:,1) + mu;
state2(:,1:3) = state2(:,1:3)*DistanceUnit;
state3(:,1) = state3(:,1) + mu;
state3(:,1:3) = state3(:,1:3)*DistanceUnit;
state4(:,1) = state4(:,1) + mu;
state4(:,1:3) = state4(:,1:3)*DistanceUnit;


for i = 1:numel(t1)
    theta = theta0 + t1(i)*theta_dot;
    T_IR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
    state1_I(i,:) = T_IR*state1(i,1:3)';
end
for i = 1:numel(t2)
    theta = theta0 + t2(i)*theta_dot;
    T_IR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
    state2_I(i,:) = T_IR*state2(i,1:3)';
end
for i = 1:numel(t3)
    theta = theta0 + t3(i)*theta_dot;
    T_IR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
    state3_I(i,:) = T_IR*state3(i,1:3)';
end
for i = 1:numel(t4)
    theta = theta0 + t4(i)*theta_dot;
    T_IR = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
    state4_I(i,:) = T_IR*state4(i,1:3)';
end


%rotating frame plots
figure
plot3(state1(:,1),state1(:,2),state1(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory a, Rotating Frame')
legend('Trajectory','Earth','Moon')

figure
plot3(state2(:,1),state2(:,2),state2(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory b, Rotating Frame')
legend('Trajectory','Earth','Moon')

figure
plot3(state3(:,1),state3(:,2),state3(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory c, Rotating Frame')
legend('Trajectory','Earth','Moon')

figure
plot3(state4(:,1),state4(:,2),state4(:,3))
hold on
plot3(-mu*DistanceUnit,0,0,'-r*')
plot3((1-mu)*DistanceUnit,0,0,'-g*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory d, Rotating Frame')
legend('Trajectory','Earth','Moon')

%inertial frame plots
figure
plot3(state1_I(:,1),state1_I(:,2),state1_I(:,3))
hold on
plot3(0,0,0,'-r*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory a, ECI Frame')
legend('Trajectory','Earth')

figure
plot3(state2_I(:,1),state2_I(:,2),state2_I(:,3))
hold on
plot3(0,0,0,'-r*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory b, ECI Frame')
legend('Trajectory','Earth')

figure
plot3(state3_I(:,1),state3_I(:,2),state3_I(:,3))
hold on
plot3(0,0,0,'-r*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory c, ECI Frame')
legend('Trajectory','Earth')

figure
plot3(state4_I(:,1),state4_I(:,2),state4_I(:,3))
hold on
plot3(0,0,0,'-r*')
grid on
grid minor
axis equal
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
title('Nondimensional Trajectory d, ECI Frame')
legend('Trajectory','Earth')
