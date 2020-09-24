clear all; close all; clc
options = [5 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 .1 0 0;0 0 0 0 .1 0;0 0 0 0 0 .1];
for i = 1:6
    global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 X Y Z K feedback
    %Constants for every scenario
    g = 9.81; %acceleration due to gravity, m/s^2
    m = .068; %kg
    rad = .06; %m
    I_x = 6.8e-5; %Components of moment of inertia matrix
    I_y = 9.2e-5;
    I_z = 1.35e-4;
    k = .0024;
    alpha = 2e-6;
    eta = 1e-3;
    K = .004;
    feedback = 0;

    f = -m*g;
    f1 = f/4;
    f2 = f/4;
    f3 = f/4;
    f4 = f/4;

    clear t y
    u_E = 0;
    v_E = 0; %m/s
    w_E = 0;
%     X = 0 + -eta*u_E^2;
%     Y = 0 + -eta*v_E^2;
%     Z = -m*g + -eta*w_E^2;
    p = options(4,i);
    q = options(5,i);
    r = options(6,i);
    phi = deg2rad(options(1,i));
    theta = deg2rad(options(2,i));
    psi = deg2rad(options(3,i));
    N = 0;
    E = 0;
    D = -1;
    vel = [u_E v_E w_E];
    omega = [p q r];
    euler = [psi theta phi];
    pos = [N E D];
    inish_condish = [vel omega euler pos 1];
    tspan = [0 5]; 
    [t,y] = ode45('odequad',tspan,inish_condish);

    figure
    subplot(1,2,1)
    hold on
    plot(t,y(:,1))
    plot(t,y(:,2))
    plot(t,y(:,3))
    title('Velocity Vs Time non-linear model, r=.1 rad/s')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    legend('u^E','v^E','w^E')


    %Assignment 3, linearized model of quadcopter dynamics and implementation
    %of control laws
    clc;
    global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4
    %Constants for every scenario
    g = 9.81; %acceleration due to gravity, m/s^2
    m = .068; %kg
    rad = .06; %m
    I_x = 6.8e-5; %Components of moment of inertia matrix
    I_y = 9.2e-5;
    I_z = 1.35e-4;
    k = .0024;
    alpha = 2e-6;
    eta = 1e-3;
    %linearized model for steady hover trim state
    f = -m*g;
    f1 = f/4;
    f2 = f/4;
    f3 = f/4;
    f4 = f/4;
    %Treat initial conditions as disturbances in "top 8" equations
    dp = options(4,i);
    dq = options(5,i);
    dr = options(6,i);
    dphi = deg2rad(options(1,i));
    dtheta = deg2rad(options(2,i));
    du = 0;
    dv = 0;
    dw = 0;
    % inish_condish = [d_Lc d_Mc d_Nc d_p d_q d_phi d_theta d_Zc];
    inish_condish = [dp dq dr dphi dtheta du dv dw];
    tspan = [0 5];
    [t y] = ode45('hw3ode',tspan,inish_condish);

    subplot(1,2,2)
    hold on
    plot(t,y(:,6))
    plot(t,y(:,7))
    plot(t,y(:,8))
    title('Velocity Vs Time linear model, r=.1 rad/s')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    legend('u^E','v^E','w^E')

end

%Derivative control law implemented
options = [.1 0 0;0 .1 0;0 0 .1];
for i = 1:3
    global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 X Y Z K feedback
    %Constants for every scenario
    g = 9.81; %acceleration due to gravity, m/s^2
    m = .068; %kg
    rad = .06; %m
    I_x = 6.8e-5; %Components of moment of inertia matrix
    I_y = 9.2e-5;
    I_z = 1.35e-4;
    k = .0024;
    alpha = 2e-6;
    eta = 1e-3;
    K = .004;
    feedback = 1;

    f = -m*g;
    f1 = f/4;
    f2 = f/4;
    f3 = f/4;
    f4 = f/4;

    clear t y
    u_E = 0;
    v_E = 0; %m/s
    w_E = 0;
%     X = 0 + -eta*u_E^2;
%     Y = 0 + -eta*v_E^2;
%     Z = -m*g + -eta*w_E^2;
    p = options(1,i);
    q = options(2,i);
    r = options(3,i);
    phi = 0;
    theta = 0;
    psi = 0;
    N = 0;
    E = 0;
    D = -1;
    vel = [u_E v_E w_E];
    omega = [p q r];
    euler = [psi theta phi];
    pos = [N E D];
    inish_condish = [vel omega euler pos 1];
    tspan = [0 .5]; 
    [t,y] = ode45('odequad',tspan,inish_condish);

    figure
    subplot(1,2,1)
    hold on
    plot(t,y(:,1))
    plot(t,y(:,2))
    plot(t,y(:,3))
    title('Velocity Vs Time, Derivative Control Law')
    xlabel('Time (s)')
    ylabel('m/s')
    legend('u^E','v_E','w_E')
    subplot(1,2,2)
    hold on
    plot(t,y(:,4))
    plot(t,y(:,5))
    plot(t,y(:,6))
    title('Yaw Rate Vs Time')
    xlabel('Time (s)')
    ylabel('rad/s')
    legend('p','q','r')
end

%comparing spidercopter reaction with and without feedback controls 

data1 = load('M_1115_A2.mat');
data2 = load('W_1134_A3.mat');

figure
hold on
plot3(data1.rt_estim.signals.values(:,1),data1.rt_estim.signals.values(:,2),-data1.rt_estim.signals.values(:,3))
plot3(data2.rt_estim.signals.values(:,1),data2.rt_estim.signals.values(:,2),-data2.rt_estim.signals.values(:,3))
title('Spiderbot Copter Position')
xlabel('N (m)')
ylabel('E (m)')
zlabel('D (m)')
legend('No Controls','Derivative Control')
figure
subplot(1,2,1)
hold on
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,10))
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,11))
plot(data1.rt_estim.time(:),data1.rt_estim.signals.values(:,12))
title('SpiderCopter No Controls')
xlabel('time (s)')
ylabel('rad/s')
legend('p','q','r')
subplot(1,2,2)
hold on
plot(data2.rt_estim.time(:),data2.rt_estim.signals.values(:,10))
plot(data2.rt_estim.time(:),data2.rt_estim.signals.values(:,11))
plot(data2.rt_estim.time(:),data2.rt_estim.signals.values(:,12))
title('SpiderCopter Derivative Control')
xlabel('time (s)')
ylabel('rad/s')
legend('p','q','r')



