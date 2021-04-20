% HW 8 -- Use single shooting algorithm to fina an orbit that is periodic
% within tolerance of 1e-12. Update z and ydot components of ICs and leave
% initial value of x unchanged. Provide plot of orbit. 
clear all;close all;clc
opts = odeset('RelTol',10e-12,'AbsTol',10e-12,'Events',@eventfunc);
% example ode call with event function
% [t,y,~,~,~] = ode45(@(t,y) odefunP6(t,y,mu),tspan,x0,optsnew);


mu = 0.012150585609624;
DistanceUnit = 384747.962856037;

ICs = [1.14 0 -0.16 0 -0.223 0]';
x = ICs;
Phi0 = eye(6);
Phi0_flat = reshape(Phi0,36,1);
tspan = [0 2*pi];
its = 1;
figure
hold on
while 1
    Z = [ICs(:,its);Phi0_flat];

    [t,x] = ode45(@(t,Z) CRTBP_ode(Z),tspan,Z,opts);
    plot(x(:,1),x(:,2))
    if  abs(x(end,4)) < 1e-8 && abs(x(end,6))
        break;
    else
        Phi = reshape(x(end,7:end),6,6);
        R1 = sqrt((x(end,1) + mu)^2 + x(end,2)^2 + x(end,3)^2);
        R2 = sqrt((x(end,1) - 1 + mu)^2 + x(end,2)^2 + x(end,3)^2);
        xdd = 2*x(end,5) + x(end,1) - (1 - mu)*((x(end,1)+mu)/(R1^3)) - mu*((x(end,1) - 1 + mu)/(R2^3));
        zdd = -(1 - mu)*(x(end,3)/(R1^3)) - mu*(x(end,3)/(R2^3));

        dxdot = -x(end,4);
        dzdot = -x(end,6);
        dX0 = inv([Phi(4,3) Phi(4,5);Phi(6,3) Phi(6,5)] - 1/x(end,5)*[xdd;zdd]*[Phi(2,3) Phi(2,5)])*[dxdot;dzdot];
        ICs(:,its+1) = ICs(:,its); 
        ICs(3,its+1) = ICs(3,its) + dX0(1);
        ICs(5,its+1) = ICs(5,its) + dX0(2);
        
        its = its + 1;
    end
end
tspan = [0 2*t(end)];
opts = odeset('RelTol',10e-12,'AbsTol',10e-12);
Z = [ICs(:,its);Phi0_flat];

[t,x] = ode45(@(t,Z) CRTBP_ode(Z),tspan,Z,opts);

figure
plot3(x(:,1),x(:,2),x(:,3),'LineWidth',1)
hold on
plot3(1-mu,0,0,'*','MarkerSize',15)
plot3(1.155682160290809,0,0,'*','MarkerSize',15)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
grid minor
legend('Trajectory','Moon','L2','FontSize',18)
set(gca,'FontSize',20)
title('Final Periodic Orbit')

figure
hold on
plot(x(:,1),x(:,2),'Linewidth',1)
hold on
plot(1-mu,0,'*','MarkerSize',15)
plot(1.155682160290809,0,'*','MarkerSize',15)
axis equal
xlabel('X')
ylabel('Y')
grid on
grid minor
legend('Trajectory','Moon','L2','FontSize',18)
set(gca,'FontSize',20)
title('Final Periodic Orbit - XY Plane')

figure
hold on
plot(x(:,1),x(:,3),'Linewidth',1)
hold on
plot(1-mu,0,'*','MarkerSize',15)
plot(1.155682160290809,0,'*','MarkerSize',15)
axis equal
xlabel('X')
ylabel('Z')
grid on
grid minor
legend('Trajectory','Moon','L2','FontSize',18)
set(gca,'FontSize',20)
title('Final Periodic Orbit - XZ Plane')


figure
hold on
plot(x(:,2),x(:,3),'Linewidth',1)
hold on
plot(0,0,'*','MarkerSize',15)
plot(0,0,'*','MarkerSize',15)
axis equal
xlabel('Y')
ylabel('Z')
grid on
grid minor
legend('Trajectory','Moon','L2','FontSize',18)
set(gca,'FontSize',20)
title('Final Periodic Orbit - YZ Plane')



