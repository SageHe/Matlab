%% Problem 3
clear all;close all;clc

tspan = [0 5];
y0 = [0 0 1000];
du = [10 100 300];
for i = 1:3
    deltau = du(i);
    [t,y] = ode45(@(t,y)odefunNL(t,y,deltau),tspan,y0);

    [t2,y2] = ode45(@(t,y)odefunL(t,y,deltau),tspan,y0);
    
    figure(i)
    hold on
    plot(t,y(:,1))
    plot(t2,y2(:,1))
    title('Rocket Altitude VS Time')
    xlabel('Time (s)')
    ylabel('Altitude (m)')
    legend('Non-linearized Solution','Linearized Solution')
end
    