clear all;close all;clc
%%Problem 1
mu = 0.1;
y0 = [.4 0 0 0];
tspan = [0 2*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
grid on
grid minor
title('Planar CR3BP trajectory for x0 = [.4 0 0 0]')

V = @(x,yv) .5.*(x.^2+yv.^2) + (1-mu)./sqrt((x+mu).^2+yv.^2) + mu./sqrt((x-1+mu).^2+yv.^2);

xprime = y(:,3);
yprime = y(:,4);

J = .5*(xprime.^2+yprime.^2) - V(y(:,1),y(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2
% clear
% x = linspace(-5,10,100000);
y = 0;
mu = [1e-5 1e-3 5e-2 1e-2 .01 .05 .1 .2 .3 .45];

% Vorig = @(x,y) .5.*(x.^2+y.^2) + (1-mu)./sqrt((x+mu).^2+y.^2) + mu./sqrt((x-1+mu).^2+y.^2);
% V = @(x) .5.*(x.^2) + (1-mu)./sqrt((x+mu).^2) + mu./sqrt((x-1+mu).^2);
for i = 1:numel(mu)
    Vx = @(x,y) x + ((mu(i) - 1).*(x+mu(i)))./(((x+mu(i)).^2 + y^2).^(3/2)) - (mu(i).*(x - 1 + mu(i)))./(((x - 1 + mu(i)).^2 + y^2).^(3/2));
    Vxs = @(x) Vx(x,0);
    soln(i,:) = fsolve(Vxs,[-1 .5 1]);
end

figure
hold on
plot(mu,soln(:,1))
plot(mu,soln(:,2))
plot(mu,soln(:,3))
grid on
grid minor
xlabel('\mu')
ylabel('X')
title('L1, L2, and L3 Points for Varying \mu')
legend('L3','L1','L2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3
% clear
mu = .5*(1-sqrt(23/27)) - 1e-2;
fun = @root2d;
L45 = fsolve(fun,[.5,.8]);

y0 = [L45(1) L45(2) 0 0];
tspan = [0 30*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
grid on
grid minor
title('Stable L_{4,5} Trajectory with x0 = [0.4715 0.8660 0 0]')

mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable Motion About L_{4,5}')

mu = .5*(1-sqrt(23/27)) + 1e-2;
fun2 = @root2d2;
L45_2 = fsolve(fun2,[.5,.8]);

y0 = [L45_2(1) L45_2(2) 0 0];
tspan = [0 20*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
grid on
grid minor
title('Unstable L_{4,5} Trajectory for x0 = [0.4515 0.8660 0 0]')

clear a e w_tilde
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Unstable Motion About L_{4,5}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 4 % 5
mu = 0.01;
L3 = soln(5,1);
L1 = soln(5,2);
L2 = soln(5,3);
y0 = [(L1+1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);
%Plot manifold
L1fig = figure;
subplot(2,2,1)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L1 Unstable Manifold, X_0 = L1+1e-8')
grid on
grid minor
%Calculate inertial pos. and vel. then orbital elements
clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end

figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Oribtal Elements for Unstable L1 point, Right Perturbed')

mu = 0.01;
y0 = [(L1-1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L1fig)
subplot(2,2,2)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L1 Unstable Manifold, X_0 = L1-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Oribtal Elements for Unstable L1 point, Left Perturbed')

y0 = [(L1+1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L1fig)
subplot(2,2,3)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L1 Stable Manifold, X_0 = L1+1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L1 point, Right Perturbed')

y0 = [(L1-1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L1fig)
subplot(2,2,4)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L1 Stable Manifold, X_0 = L1-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L1 point, Left Perturbed')
%L2 point
y0 = [(L2+1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

L2fig = figure;
subplot(2,2,1)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L2 Unstable Manifold, X_0 = L2+1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Unstable L2 point, Right Perturbed')

y0 = [(L2-1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L2fig)
subplot(2,2,2)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L2 Unstable Manifold, X_0 = L2-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Unstable L2 point, Left Perturbed')

y0 = [(L2+1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L2fig)
subplot(2,2,3)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L2 Stable Manifold, X_0 = L2+1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L2 point, Right Perturbed')

y0 = [(L2-1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L2fig)
subplot(2,2,4)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L2 Stable Manifold, X_0 = L2-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L2 point, Left Perturbed')
%L3 points
y0 = [(L3+1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

L3fig = figure;
subplot(2,2,1)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L3 Unstable Manifold, X_0 = L3+1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Unstable L3 point, Right Perturbed')

y0 = [(L3-1e-8) 0 0 0];
tspan = [0 4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L3fig)
subplot(2,2,2)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L3 Unstable Manifold, X_0 = L3-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Unstable L3 point, Left Perturbed')

y0 = [(L3+1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L3fig)
subplot(2,2,3)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L3 Stable Manifold, X_0 = L3+1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L3 point, Right Perturbed')

y0 = [(L3-1e-8) 0 0 0];
tspan = [0 -4*pi];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) odefun(t,y,mu),tspan,y0,opts);

figure(L3fig)
subplot(2,2,4)
hold on
plot(y(:,1),y(:,2))
xlabel('X')
ylabel('Y')
title('L3 Stable Manifold, X_0 = L3-1e-8')
grid on
grid minor

clear a e w_tilde pos_I vel_I T
mu = 1;
for i = 1:numel(t)
    T = [cos(t(i)) -sin(t(i));sin(t(i)) cos(t(i))];
    pos_I = T*y(i,[1:2])';
    pos_I = [pos_I;0];
    vel_I = T*[(y(i,3)-y(i,2)) ; (y(i,4)+y(i,1))];
    vel_I = [vel_I;0];
    [a(i),e(i),w_tilde(i),~,~,~,~,~,~] = compOE(pos_I,vel_I,mu);
end
mu = 0.01;
figure
subplot(3,1,1)
hold on
plot(t,a)
xlabel('Time (s)')
ylabel('Semi-Major Axis (a)')
grid on
grid minor
subplot(3,1,2)
plot(t,e)
xlabel('Time (s)')
ylabel('Eccentricity (e)')
grid on
grid minor
subplot(3,1,3)
plot(t,w_tilde)
xlabel('Time (s)')
ylabel(' ($\tilde w$)','Interpreter','Latex')
grid on 
grid minor
sgtitle('Orbital Elements for Stable L3 point, Left Perturbed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 6
mu = 0.01;
optsnew = odeset('RelTol',1e-10,'AbsTol',1e-12,'Events',@eventfunc);
stm0 = reshape(eye(4),16,1);
x0 = [L1-2.51e-2;0;0;0.2510;stm0];
% y0 = [(L1+1e-8) 0 0 0];
tspan = [0 2*pi];
% opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y,~,~,~] = ode45(@(t,y) odefunP6(t,y,mu),tspan,x0,optsnew);
% figure
% hold on
% plot(y(:,1),y(:,2))
err_tol = 1e-12;
iter_max = 100;
yd0_new = x0(4);
err_dxf = 100;
iter = -1;
figure
hold on
while (err_dxf > err_tol) && (iter <= iter_max)
    iter = iter + 1;
    
    if iter > 0
        gamma = stm_tf_t0(1,4);
        
       dyd0 = (gamma^-1)*dxf;
       
       yd0_new = x0_new(4) - dyd0;
    end
    
    x0_new = [x0(1:3);yd0_new;stm0];
    
    [t,x,~,~,~] = ode45(@(t,y) odefunP6(t,y,mu),tspan,x0_new,optsnew);
%     plot(x(:,1),x(:,2))
    
    stm_tf_t0 = reshape(x(end,5:20),4,4);

    dxf = x(end,1) - x0(1);
    err_dxf = abs(dxf);
end
plot(x(:,1),x(:,2))
xlabel('X')
ylabel('Y')
title('Lyapunov PO about L1, Propogated Backwards')
grid on
grid minor
fprintf('%1d iterations ... Error: %1.2e\n', iter, err_dxf)
% plot(x(:,1),x(:,2))
