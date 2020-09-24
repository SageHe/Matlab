%% Simulate nonlinear dynamical system to obtain nominal trajectory
clear;close all;clc

mu = 398600;
r = 6678;

X0 = [6678 0 0 r*sqrt(mu/r^3)]';
tvec = linspace(0,14000,1401);
tspan = [0 14000];
opts = odeset('RelTol',1e-13);
[t,X] = ode45('odeFunc',tspan,X0,opts);
X = interp1(t,X,tvec);
X = X';
%simulate nonlinear model and use perturbations for ground truth model
X0 = [6678 0 0 r*sqrt(mu/r^3)]';
pert = [0 0.075 0 -0.021]';
X0 = X0+pert;
tvec = linspace(0,14000,1401);
tspan = [0 14000];
opts = odeset('RelTol',1e-13);
[t,X_gt] = ode45('odeFunc',tspan,X0,opts);
X_gt = interp1(t,X_gt,tvec);
X_gt = X_gt';

% mu = 398600;
% r = 6678;
% A = [0 1 0 0;-mu/(R^3) 0 0 0;0 0 0 1;0 0 -mu/(R^3) 0];
B = [0 0;1 0;0 0;0 1];
L = B;
dt = 10; %seconds
% 
% Ahat = [A B;zeros(2,6)];
% eAhat = expm(Ahat*dt);
% F = eAhat([1:4],[1:4]);
% G = eAhat([1:4],[5:6]);
% 
% [V,D] = eig(F);
% 
% P = [G F*G F^2*G F^3*G];
% 
% rank(P)
%% part c 
% Simulate linearized DT dynamics and measurement models near lin. point
% near for system, assuming a reasonable initial state pert. from the
% linearization point and assuming no process noise, measurement noise, or
% control input perts. Use results to compare and validate your Jacobians
% and DT against a full nonlin. sim. of the system dynamics and
% measurements using ode45 and provide plots to compare linearized DT and
% nonlin. DT model. 
%Define given constants
Re = 6378;
We = (2*pi)/86400;
X0 = [6678 0 0 (r)*sqrt(mu/(r^3))]';
pert = [0 0.075 0 -0.021]';
Xlin = [X0+pert];
%Determine H matrix and calculate observations y over one orbit period
y = [];
Y = [];
for i = 1:1400
    Ytemp = [];
    for j = 1:12
        Ytemp = [Ytemp;nom_meas(X(:,i),X_gt(:,i),tvec(i),j)];
    end
    Y = [Y Ytemp];
end
F = [];
y = Y(:,1);
y([4:end],1) = nan;
for i = 1:1400
    F(:,:,i) = fFromA(X(1,i),X(3,i));
    Xlin(:,i+1) = X(:,i+1) + F(:,:,i)*(Xlin(:,i)-X(:,i));
    y_loop = [];
    Yloop = [];
    Ytemp = [];
    for j = 1:12
        theta0(j) = (j - 1)*(pi/6);
        Xs = Re*cos(We*(i*10) + theta0(j));
        Xdots = -Re*sin(We*(i*10)+theta0(j))*We;
        Ys = Re*sin(We*(i*10) + theta0(j));
        Ydots = Re*cos(We*(i*10) + theta0(j))*We; 
        rho = sqrt((X(1,i)-Xs)^2 + (X(3,i)-Ys)^2);
        rhodot = ((X(1,i) - Xs)*(X(2,i)-Xdots) + (X(3,i)-Ys)*(X(4,i)-Ydots))/rho;
        num = ((X(1,i) - Xs)*(X(2,i)-Xdots) + (X(3,i)-Ys)*(X(4,i)-Ydots));
        phi = atan2((X(3,i)-Ys),(X(1,i)-Xs));
        phi_check = atan2((X_gt(3,i)-Ys),(X_gt(1,i)-Xs));
        theta_i = atan2(Ys,Xs);
        %calculate nominal measurements
%         Ytemp = [Ytemp; nom_meas(X(:,i),tvec(i),j)];
        H(:,:,i,j) = [(X(1,i)-Xs)/rho 0 (X(3,i)-Ys)/rho 0;
                (rho*(X(2,i)-Xdots)-(num)*((X(1,i)-Xs)/rho))/(rho^2) (X(1,i)-Xs)/rho  (rho*(X(4,i)-Ydots)-(num)*((X(3,i)-Ys)/rho))/(rho^2) (X(3,i)-Ys)/rho;
                (Ys-X(3,i))/(rho^2) 0 (X(1,i)-Xs)/(rho^2) 0];
        lower1 = -pi/2+theta_i;
        upper1 = pi/2 +theta_i;
        lower2 = lower1;
        upper2 = upper1;
        if -pi/2+theta_i < -pi
            lower1 = -pi/2+theta_i+2*pi;
            upper1 = pi;
            lower2 = -pi;
            upper2 = pi/2+theta_i;
        else
        if  pi/2+theta_i > pi
            lower1 = -pi;
            upper1 = pi/2+theta_i-2*pi;
            lower2 = -pi/2+theta_i;
            upper2 = pi;
        end
        end
        if (phi_check >= lower1 && phi_check <= upper1)|| (phi_check >= lower2 && phi_check <= upper2)
            y_temp = Y([3*j-2:3*j],i) + H(:,:,i,j)*(Xlin(:,i)-X(:,i));
        else
            y_temp = NaN(3,1);
        end
%         if phi >= ((-pi/2)+atan2(Ys,Xs)) && phi <= ((pi/2)+atan2(Ys,Xs))
%             ytemp = H*(X(:,i) - [Xs Xdots Ys Ydots]');
%         else
%             ytemp = nan(3,1);
%         end
        y_loop = [y_loop; y_temp];
%         Yloop = [Yloop; Ytemp];
    end
    y = [y y_loop];
%     Y = [Y Ytemp];
end
figure
subplot(2,2,1)
% plot(tvec,Xlin(1,:))
% hold on
plot(tvec,X(1,:))
subplot(2,2,2)
% plot(tvec,Xlin(2,:))
% hold on
plot(tvec,X(2,:))
subplot(2,2,3)
% plot(tvec,Xlin(3,:))
% hold on
plot(tvec,X(3,:))
subplot(2,2,4)
% plot(tvec,Xlin(4,:))
% hold on
plot(tvec,X(4,:))
figure
subplot(2,2,1)
plot(tvec,Xlin(1,:) - X(1,:))
title('X')
subplot(2,2,2)
plot(tvec,Xlin(2,:) - X(2,:))
title('Xdot')
subplot(2,2,3)
plot(tvec,Xlin(3,:) - X(3,:))
title('Y')
subplot(2,2,4)
plot(tvec,Xlin(4,:) - X(4,:))
title('Ydot')
% tvec(end) = [];
figure
subplot(3,1,1)
for i=1:12
    scatter(tvec,y(3*i-2,:))
    hold on
end
subplot(3,1,2)
for i=1:12
    scatter(tvec,y(3*i-1,:))
    hold on
end
subplot(3,1,3)
for i=1:12
    scatter(tvec,y(3*i,:))
    hold on
end
D = load('orbitdeterm_finalproj_KFdata.mat');
Omega = dt*L;
e = (X(:,2) - Xlin(:,2))*(X(:,2) - Xlin(:,2))';
e = diag(e);
e = diag(e);
Pk_p = cov(e);
Pk_p(2,2) = Pk_p(2,2)*.05;
D.Qtrue = D.Qtrue*.01;
D.Qtrue(1,1) = D.Qtrue(1,1)*0.1;
D.ydata(1) = [];
D.tvec(1) = [];
Sv = chol(D.Rtrue,'lower');
rng(100)
% LKF loop
for i = 1:numel(D.ydata)
    if ~isempty(D.ydata{i})
        D.ydata{i}(1:3,:) = D.ydata{i}(1:3,:) + Sv*randn(3,1);
    end
end
for i = 1:1400
    dx_kp1_m = Xlin(:,i)-X(:,i);
%     P_kp1_m = F(:,:,i)*Pk_p*F(:,:,i)' + Omega*D.Qtrue*Omega';
    if isempty(D.ydata{i})
    dy_kp1(:,i) = dy_kp1(:,i-1);
    dx_kp1_p(:,i) = dx_kp1_p(:,i-1);
    P_kp1_p(:,:,i) = P_kp1_p(:,:,i-1);
    continue 
    end
    stat_num = DST(Y(:,i),D.ydata{i}(end,:));
    for j = 1:size(stat_num)
        s = stat_num(j);
        P_kp1_m = F(:,:,i)*Pk_p*F(:,:,i)' + Omega*D.Qtrue*Omega';
        K_kp1 = P_kp1_m*H(:,:,i,s)'*inv(H(:,:,i,s)*P_kp1_m*H(:,:,i,s)' + D.Rtrue);
        dy_kp1(:,i) = D.ydata{i}(1:3,j) - Y([3*s-2:3*s],i);
        dx_kp1_p(:,i) = dx_kp1_m + K_kp1*(dy_kp1(:,j) - H(:,:,i,s)*dx_kp1_m);
        P_kp1_p(:,:,i) = (eye(4) - K_kp1*H(:,:,i,s))*P_kp1_m;
        Pk_p = P_kp1_p(:,:,i);
        if numel(stat_num) > 1
            dx_kp1_m = dx_kp1_p;
        end
    end
end
X_KF = X(:,[1:1400]) + dx_kp1_p;
Sv = chol(D.Qtrue*1e7,'lower');
X_gt = X_gt + Omega*Sv*randn(2,1);
figure
subplot(4,1,1)
plot(D.tvec,dx_kp1_p(1,:))
subplot(4,1,2)
plot(D.tvec,dx_kp1_p(2,:))
subplot(4,1,3)
plot(D.tvec,dx_kp1_p(3,:))
subplot(4,1,4)
plot(D.tvec,dx_kp1_p(4,:))
X_gt(:,end) = [];
figure
sgtitle('States of Spacecraft')
subplot(2,2,1)
plot(D.tvec, X_KF(1,:))
hold on
plot(D.tvec, X_gt(1,:))
ylabel('X Position (km)')
xlabel('Time (s)')
title('$X$','interpreter','latex')
subplot(2,2,2)
plot(D.tvec, X_KF(2,:))
hold on
plot(D.tvec, X_gt(2,:))
legend('KF','Ground Truth')
ylabel('X Velocity (km/s)')
xlabel('Time (s)')
title('$\dot{X}$','interpreter','latex')
subplot(2,2,3)
plot(D.tvec, X_KF(3,:))
hold on
plot(D.tvec, X_gt(3,:))
ylabel('Y Position (km)')
xlabel('Time (s)')
title('$Y$','interpreter','latex')
subplot(2,2,4)
plot(D.tvec, X_KF(4,:))
hold on
plot(D.tvec, X_gt(4,:))
ylabel('Y Velocity (km/s)')
xlabel('Time (s)')
title('$\dot{Y}$','interpreter','latex')
figure
subplot(2,2,1)
plot(D.tvec,X_gt(1,:) - X_KF(1,:))
subplot(2,2,2)
plot(D.tvec,X_gt(2,:) - X_KF(2,:))
subplot(2,2,3)
plot(D.tvec,X_gt(3,:) - X_KF(3,:))
subplot(2,2,4)
plot(D.tvec,X_gt(4,:) - X_KF(4,:))
% EKF loop
e = (X(:,2) - Xlin(:,2))*(X(:,2) - Xlin(:,2))';
e = diag(e);
e = diag(e);
Pk_p = cov(e);
Pk_p(4,4) = Pk_p(4,4)*.05;
D.Qtrue = D.Qtrue*.01;
D.Qtrue(1,1) = D.Qtrue(1,1)*0.1;
Xk_p = X(1,:); %should this be X_gt?
X0 = [6678 0 0 r*sqrt(mu/r^3)]';
pert = [0 0.075 0 -0.021]';
X0 = X0+pert;
tspan = [0 10];
% y_kp1_m = Y([1:3],1);
for i = 1:1400
    [t,X_GT] = ode45('odeFunc',tspan,X0,opts);
    x_kp1_m = X_GT(end,:)';
%     P_kp1_m = F(:,:,i)*Pk_p*F(:,:,i)' + Omega*D.Qtrue*Omega';
    if isempty(D.ydata{i})
    y_kp1_m(:,i) = y_kp1_m(:,i-1);
    x_kp1_p(:,i) = x_kp1_p(:,i-1);
    P_kp1_p(:,:,i) = P_kp1_p(:,:,i-1);
    continue 
    end
    stat_num = DST(Y(:,i),D.ydata{i}(end,:)); %Pretty sure this is using the wrong nonlin mesasurements currently 
    for j = 1:size(stat_num,2)
        s = stat_num(j);
        F(:,:,i) = fFromA(x_kp1_m(1,1),x_kp1_m(3,1));
        P_kp1_m = F(:,:,i)*Pk_p*F(:,:,i)' + Omega*D.Qtrue*Omega';
        H(:,:,i,s) = calcH(x_kp1_m,i,s);
        K_kp1 = P_kp1_m*H(:,:,i,s)'*inv(H(:,:,i,s)*P_kp1_m*H(:,:,i,s)' + D.Rtrue);
        y_kp1_m(:,i) = nom_meas_EKF(x_kp1_m,X_GT(end,:),tvec(i),s);
        e_y_kp1 = D.ydata{i}([1:3],j) - y_kp1_m(:,i);
        x_kp1_p(:,i) = x_kp1_m + K_kp1*e_y_kp1;
        P_kp1_p(:,:,i) = (eye(4) - K_kp1*H(:,:,i,s))*P_kp1_m;
        Pk_p = P_kp1_p(:,:,i);
        if numel(stat_num) > 1
            x_kp1_m = x_kp1_p(:,i);
        end
    end
    tspan = tspan + 10;
    X0 = x_kp1_p(:,i);
end
X_EKF = x_kp1_p;
figure
sgtitle('States of Spacecraft')
subplot(2,2,1)
plot(D.tvec, X_EKF(1,:))
hold on
plot(D.tvec, X_gt(1,:))
ylabel('X Position (km)')
xlabel('Time (s)')
title('$X$','interpreter','latex')
subplot(2,2,2)
plot(D.tvec, X_EKF(2,:))
hold on
plot(D.tvec, X_gt(2,:))
legend('KF','Ground Truth')
ylabel('X Velocity (km/s)')
xlabel('Time (s)')
title('$\dot{X}$','interpreter','latex')
subplot(2,2,3)
plot(D.tvec, X_EKF(3,:))
hold on
plot(D.tvec, X_gt(3,:))
ylabel('Y Position (km)')
xlabel('Time (s)')
title('$Y$','interpreter','latex')
subplot(2,2,4)
plot(D.tvec, X_EKF(4,:))
hold on
plot(D.tvec, X_gt(4,:))
ylabel('Y Velocity (km/s)')
xlabel('Time (s)')
title('$\dot{Y}$','interpreter','latex')

figure
subplot(2,2,1)
plot(D.tvec,X_gt(1,:) - X_EKF(1,:))
subplot(2,2,2)
plot(D.tvec,X_gt(2,:) - X_EKF(2,:))
subplot(2,2,3)
plot(D.tvec,X_gt(3,:) - X_EKF(3,:))
subplot(2,2,4)
plot(D.tvec,X_gt(4,:) - X_EKF(4,:))

% for i = 1:1400
%     x_kp1_m = Xk_p;
%     F(:,:,i) = fFromA(x_kp1_m(1,:),x_kp1_m(3,:));
%     P_kp1_m = F(:,:,i)*Pk_p*F(:,:,i) + Omega*10*D.Qtrue*Omega'; % may want to change alteration of D.Qtrue in LKF to separate variable so as not to mess up EKF by calling D.Qtrue again
%     for j = 1:size(D.ydata{i},2)
%         H(:,:,i,j) = [(x_kp1_m(1,i)-Xs)/rho 0 (x_kp1_m(3,i)-Ys)/rho 0;
%         (rho*(X(2,i)-Xdots)-(num)*((X(1,i)-Xs)/rho))/(rho^2) (X(1,i)-Xs)/rho  (rho*(X(4,i)-Ydots)-(num)*((X(3,i)-Ys)/rho))/(rho^2) (X(3,i)-Ys)/rho;
%         (Ys-X(3,i))/(rho^2) 0 (X(1,i)-Xs)/(rho^2) 0];
% 
%         s = D.ydata{i}(end,j);
%         y_kp1_m = Y(:,i);
%         K_kp1 = P_kp1_m*H(:,:,i,s)'*inv(H(:,:,i,s)*P_kp1_m*H(:,:,i,s)' + D.Rtrue);
%         P_kp1_p(:,:,j) = (eye(4) - K_kp1*H(:,:,i,s))*P_kp1_m;
%     end
%     if isempty(D.ydata{i})
%         dy_kp1(:,i) = dy_kp1(:,i-1);
%         dx_kp1_p_KF(:,i) = dx_kp1_p_KF(:,i-1);
%         P_kp1_p(:,:,i) = P_kp1_p(:,:,i-1);
%     else
%         dy_kp1(:,i) = mean(dy_kp1_loop,2);
%         dx_kp1_p_KF(:,i) = mean(dx_kp1_p,2);
%         P_kp1_p(:,:,i) = mean(P_kp1_p,3);
%     end
%     Pk_p = P_kp1_p(:,:,i);
% end

% Last left off at line 295, seems like y_kp1_m is being calculated
% incorrectly, resulting in nan y_kp1_m vals for the rest of the sim,
% leading to singular near working precision P matrices. Check why y_kp1_m
% is coming out as nan after iteration 28, what is skewing y_kp1_m
% calculation
    
    