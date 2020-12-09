clear all;close all;clc 
%% Problem 1 part C
% part i
% Define constants
 J2 = 1e-3;
 mu = 4e5;
 Re = 6400;
 rho = (2*pi)/(365*24*3600);
 temp = ((2*rho)/(3*J2*Re^2*sqrt(mu)))^(2/7);
 a_M = 1/temp;
 % pot curve of relation showing inc. for SSO as a func. of orbit radius
 a = linspace(6400,a_M,100);
 i = acos(-(a./a_M).^(7/2));
 
 figure
 plot(a,i)
 xlabel('Orbit Radius (Km)')
 ylabel('Inclination (rad)')
 title('Orbit Radius VS Inclination')
 grid on
 grid minor
 
 %%Problem 2 part A
a = 7000;
e = 0.1;
i = 45; %degrees
w = 0;
Omega = 0;
sigma = 0;
nhat_Omega = cosd(Omega)*[1 0 0] + sind(Omega)*[0 1 0];
nhat_Omega_perp = -cosd(i)*sind(Omega)*[1 0 0] + cosd(i)*cosd(Omega)*[0 1 0] + sind(i)*[0 0 1];
ehat = cosd(w)*nhat_Omega + sind(w)*nhat_Omega_perp;
ehat_perp = -sind(w)*nhat_Omega + cosd(w)*nhat_Omega_perp;
n = sqrt(mu/a^3);

P = a*(1 - e^2);
E = 0;
[r0,v0] = calcRV(E,P,e,ehat,ehat_perp,mu);
T = 2*pi*sqrt(a^3/mu);
params = [J2 Re mu];
tspan = [0 100*T];
y0 = [r0 v0]';
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,params),tspan,y0,opts);
% Plot orbital position over the course of 100 oribts
figure
plot3(y(:,1),y(:,2),y(:,3))
grid on
grid minor
title('Orbit Trajectory ODE')
xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('Z Position (m)')

% Calc. and plot angular momentum
zhat = [0 0 1];
h = cross(y(:,1:3),y(:,4:6));
h = h(:,3);
figure
plot(t,h)
xlabel('Time (s)')
ylabel('Angular Momentum (Km^2/s)')
title('Angular Momentum VS Time')
ylim([3.72282e4 3.723e4])
grid on
grid minor

% Calc. and plot energy
for k = 1:numel(h)
    r = norm(y(k,1:3));
    v = norm(y(k,4:6));
    z = y(k,3);
    E(k) = .5*v^2 - (mu/r) - (mu/(2*r^3))*Re^2*J2*(1 - 3*(z^2/r^2)); 
end

figure
plot(t,E)
xlabel('Time (s)')
ylabel('Energy (km^2/s^2)')
title('Energy VS Time')
ylim([-29 -28])
grid on
grid minor
% Part b

a_avg = a;
e_avg = e;
i_avg = i;

tspan = [0 10*T]; % 10 orbits
[t,y] = ode45(@(t,y) odefun(t,y,params),tspan,y0,opts);
clear k
for k = 1:size(y,1)
    [a(k),e_temp,i(k),w(k),Omega(k),~,~,~] = compOE(y(k,1:3),y(k,4:6),mu);
    e(k) = norm(e_temp);
end
% calc orbit elements from averaging theory to compare to numberically
% integrated values
a_avg = a_avg*ones(size(t,1),1);
e_avg = e_avg*ones(size(t,1),1);
i_avg = i_avg*ones(size(t,1),1);
i_avg = deg2rad(i_avg);
p = a_avg(1)*(1 - e_avg(1)^2);
nbar = n*(1 + (3/2)*((J2*Re^2)/p^2)*(1 - (3/2)*sind(i_avg(1))^2)*(1 - e_avg(1)^2)^(.5));
w_avg = ((3/2)*((J2*Re^2)/p^2)*nbar*(2 - (5/2)*sin(i_avg(1))^2)).*t;
Omega_avg = -(3/2)*((J2*Re^2)/p^2)*nbar*cos(i_avg(1)).*t;
% Plot values to compare
figure
subplot(2,3,1)
plot(t,a)
hold on
plot(t,a_avg)
ylim([a_avg(1)-10 a_avg(1)+5])
xlabel('Time (s)')
ylabel('(Km)')
title('Semi-Major Axis')
grid on
grid minor

subplot(2,3,2)
plot(t,e)
hold on
plot(t,e_avg)
ylim([e_avg(1)-.01 e_avg(1)+.01])
xlabel('Time (s)')
title('Eccentricity')
grid on
grid minor

subplot(2,3,3)
plot(t,i)
hold on
plot(t,i_avg)
ylim([i_avg(1)-.005 i_avg(1)+.005])
xlabel('Time (s)')
ylabel('rads')
title('Inclination')
grid on
grid minor

subplot(2,3,4)
plot(t,w)
hold on
plot(t,w_avg)
xlabel('Time (s)')
ylabel('rad')
title('Argument of Periapsis')
grid on
grid minor

subplot(2,3,5)
plot(t,Omega)
hold on
plot(t,Omega_avg)
xlabel('Time (s)')
ylabel('rad')
title('Right Ascension of Ascending Node')
grid on
grid minor
%% Problem 4
mu = 1;
a = 1;
i = 0;
w = 0;
e = 0;
Omega = 0;
nhat_Omega = cosd(Omega)*[1 0 0] + sind(Omega)*[0 1 0];
nhat_Omega_perp = -cosd(i)*sind(Omega)*[1 0 0] + cosd(i)*cosd(Omega)*[0 1 0] + sind(i)*[0 0 1];
ehat = cosd(w)*nhat_Omega + sind(w)*nhat_Omega_perp;
ehat_perp = -sind(w)*nhat_Omega + cosd(w)*nhat_Omega_perp;
n = sqrt(mu/a^3);

P = a*(1 - e^2);
E = 0;
[r0,v0] = calcRV(E,P,e,ehat,ehat_perp,mu);
y0 = [r0 v0];
T = 2*pi*sqrt(a^3/mu);
tspan = [0 30*T];
g = [0.001 0.01 0.1 1];
for k = 1:4
    U = [g(k) 0 0]';
    params = {mu Re U};
    [t4{k},y4{k}] = ode45(@(t,y) odefunp4(t,y,params),tspan,y0,opts);
end
% calc angular momentum and energy to show that intergrals are constant
h1 = cross(y4{1}(:,1:3),y4{1}(:,4:6));
h1 = h1(:,1);
figure
plot(t4{1},h1)
xlabel('Time (s)')
ylabel('Angular Momentum (km^2/s^2)')
title('Angular Momentum VS Time for g = 0.001')
grid on
grid minor

h2 = cross(y4{2}(:,1:3),y4{2}(:,4:6));
h2 = h2(:,1);
figure
plot(t4{2},h2)
xlabel('Time (s)')
ylabel('Angular Momentum (km^2/s^2)')
title('Angular Momentum VS Time for g = 0.01')
grid on
grid minor

h3 = cross(y4{3}(:,1:3),y4{3}(:,4:6));
h3 = h3(:,1);
figure
plot(t4{3},h3)
xlabel('Time (s)')
ylabel('Angular Momentum (km^2/s^2)')
title('Angular Momentum VS Time for g = 0.1')
grid on
grid minor

h4 = cross(y4{4}(:,1:3),y4{4}(:,4:6));
h4 = h4(:,1);
figure
plot(t4{4},h4)
xlabel('Time (s)')
ylabel('Angular Momentum (km^2/s^2)')
title('Angular Momentum VS Time for g = 1')
grid on
grid minor

%calc energy to show constant integral
for i = 1:length(t4{1})
    v = norm(y4{1}(i,4:6));
    r = norm(y4{1}(i,1:3));
    x = y4{1}(i,1);
    E(i) = .5*v^2 - (1/r) - g(1)*x;
end
figure
plot(t4{1},E)
ylim([-.502 -.5])
xlabel('Time (s)')
ylabel('Energy (km^2/s^2)')
title('Energy VS Time for g = 0.001')
grid on
grid minor
clear v r x E

for i = 1:length(t4{2})
    v = norm(y4{2}(i,4:6));
    r = norm(y4{2}(i,1:3));
    x = y4{2}(i,1);
    E(i) = .5*v^2 - (1/r) - g(2)*x;
end
figure
plot(t4{2},E)
ylim([-.52 -.5])
xlabel('Time (s)')
ylabel('Energy (km^2/s^2)')
title('Energy VS Time for g = 0.01')
grid on
grid minor
clear v r x E

for i = 1:length(t4{3})
    v = norm(y4{3}(i,4:6));
    r = norm(y4{3}(i,1:3));
    x = y4{3}(i,1);
    E(i) = .5*v^2 - (1/r) - g(3)*x;
end
figure
plot(t4{3},E)
ylim([-.7 -.5])
xlabel('Time (s)')
ylabel('Energy (km^2/s^2)')
title('Energy VS Time for g = 0.1')
grid on
grid minor
clear v r x E

for i = 1:length(t4{4})
    v = norm(y4{4}(i,4:6));
    r = norm(y4{4}(i,1:3));
    x = y4{4}(i,1);
    E(i) = .5*v^2 - (1/r) - g(4)*x;
end
figure
plot(t4{4},E)
ylim([-1.51 -1.49])
xlabel('Time (s)')
ylabel('Energy (km^2/s^2)')
title('Energy VS Time for g = 1')
grid on
grid minor

%plot trajectories
figure
plot3(y4{1}(:,1),y4{1}(:,2),y4{1}(:,3))
grid on
grid minor
xlabel('X Position (Km)')
ylabel('Y Position (Km)')
zlabel('Z Position (Km)')
title('Orbit Trajector for 30 Orbits with g = 0.001')

figure
plot3(y4{2}(:,1),y4{2}(:,2),y4{2}(:,3))
grid on
grid minor
xlabel('X Position (Km)')
ylabel('Y Position (Km)')
zlabel('Z Position (Km)')
title('Orbit Trajector for 30 Orbits with g = 0.01')

figure
plot3(y4{3}(:,1),y4{3}(:,2),y4{3}(:,3))
grid on
grid minor
xlabel('X Position (Km)')
ylabel('Y Position (Km)')
zlabel('Z Position (Km)')
title('Orbit Trajector for 30 Orbits with g = 0.1')

figure
plot3(y4{4}(:,1),y4{4}(:,2),y4{4}(:,3))
grid on
grid minor
xlabel('X Position (Km)')
ylabel('Y Position (Km)')
zlabel('Z Position (Km)')
title('Orbit Trajector for 30 Orbits with g = 1')