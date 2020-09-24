%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodore Trozinski
% Maya West
% ASEN 5044 Estimation; Final Project
% Created: April 26, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
%% Housekeeping
clear, close all
clc
%%%%%%%
% Load Provided File
load('orbitdeterm_finalproj_KFdata.mat')

%% Setup
N = 1; % Number of Orbits
dt = 10; % timestep
% Inish Condish
mu = 398600;
Re = 6378;
Xnom = 6678;
Ynom = 0;
r0 = sqrt(Xnom^2+Ynom^2);
dXnom = 0;
dYnom = r0*sqrt(mu/r0^3);
% Nominal Conditions
X = Xnom; Y = Ynom; dX = dXnom; dY = dYnom;
% Time-invariant Matrix setup
B = [0 0; 1 0; 0 0; 0 1];
Gamma = [0 0; 1 0; 0 0; 0 1];
Omega = dt*Gamma;

%% ODE for Non-Linear Dynamics NOMINAL CASE
T = 2*pi*sqrt(r0^3/mu); % Period, [s]
xvec = [Xnom dXnom Ynom dYnom]'; % X dX Y dY
% Call
t = 0;
tvec = t:dt:N*T; % time vector for nominal (NonLinear) simulation
opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
[tstar,xstar] = ode45(@(t,x) OrbitNL(t,x,mu),tvec,xvec,opts);
xstar = xstar'; tstar = tstar'; % correct orientation, now many column vectors

%% State Propagation and Observation Calculation
% 1 Setup Vectors %
n = length(xstar);
% state
x = zeros(4,n);
% perturbation
dx = zeros(4,n);
dx(:,1) = [0.01 0.001 0.01 0.001]'*1;
% Cholesky
R = chol(Rtrue);
% Measurement Stuff
% measurement_i = {}; idk how to initialize when they do change every time
% Calculate k = zero stuff:
x(:,1) = xstar(:,1) + dx(:,1); % initial state

% Simulate with Purturbation, Nonlinear truth
opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
[~,xNL] = ode45(@(t,x) OrbitNL(t,x,mu),tvec,x(:,1),opts);
xNL = xNL';

Q = chol(Qtrue);
% Process Noise Truth
xNLnoise(:,1) = xNL(:,1);
for k = 1:(length(xstar)-1)
    opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
    noise = Q*randn(2,1);
    tvec = [10*(k-1) 10*k];
    [~,xNLnoisek] = ode45(@(t,x) OrbitNLnoise(t,x,mu,noise),tvec,xNLnoise(:,k),opts);
    xNLnoise(:,k+1) = xNLnoisek(end,:)';
end

% For loop for State Propegation and Observation Calculation
for k = 1:(n-1) % LINEARIZED DYNAMICS
    % Solve for A
    X = xstar(1,k);
    Y = xstar(3,k);
    A = [0 1 0 0;
        mu*(2*X^2-Y^2)/((Y^2+X^2).^(5/2)) 0 3*mu*X*Y/((X^2+Y^2).^(5/2)) 0;
        0 0 0 1;
        3*mu*X*Y/((X^2+Y^2).^(5/2)) 0 mu*(2*Y^2-X^2)/((Y^2+X^2).^(5/2)) 0];
    F{k} = eye(4) + dt*A;
    
    % propegate perturbation
    dx(:,k+1) = F{k}*dx(:,k);% + Omega*w(k);
    % Linear State Calculation
    x(:,k+1) = xstar(:,k+1) + dx(:,k+1);
    % Measurements
    t = (k-1)*10;
    [H{k+1},included{k+1},ystar_{k+1},ystar{k+1}] = assembleH(t,xstar(:,k+1));
    
    % Nonlinear recieved measurements
    [~,~,~,ytruth{k+1}] = assembleH(t,xNLnoise(:,k+1));
    
    % Actual y:
    %[~,~,y{k+1},~] = assembleH(t,x(:,k+1));
    yNL{k+1} = ytruth{k+1};
    % Add Observation Noise
    rng(100)
    if numel(ytruth{k+1}) == 4
        yNL{k+1}(1:3) = ytruth{k+1}(1:3) + R*randn(3,1);
    end
    if numel(ytruth{k+1}) > 4
        yNL{k+1}(1:3,:) = ytruth{k+1}(1:3,:) + R*randn(3,2);
    end
end


%% Plots:
% Plot States:
figure(1) % plot of nominal States, X
subplot(2,1,1)
plot(tstar,xNL(1,:))
hold on
plot(tstar,x(1,:))
legend('Non-Linear','Linearized')
title('X Altitude over time')
xlabel('Time, [s]')
ylabel('Distance, [km]')
subplot(2,1,2)
plot(tstar,xNL(2,:))
hold on
plot(tstar,x(2,:))
legend('Non-Linear','Linearized')
title('X Velocity over time')
xlabel('Time, [s]')
ylabel('Velocity, [km/s]')

figure(2) % plot of error from linearized to NL simulation
subplot(2,1,1)
plot(tstar,x(1,:)-xNL(1,:))
hold on
title('X Positional Error Between Linearized and Nonlinear Simulations')
xlabel('Time, [s]')
ylabel('Distance Error, [km]')
subplot(2,1,2)
plot(tstar,x(2,:)-xNL(2,:))
hold on
title(' X Velocity Error Between Linearized and Nonlinear Simulations')
xlabel('Time, [s]')
ylabel('Velocity Error, [km/s]')

%% NEES NIS Test Setup
% Set up Monte Carlo
rng(100)
Nruns = 1; % number of monte carlo sim runs <- takes a long time with many MCS
xvars = 4; % not sure what this is, guessing these are num of variables in H, one for 
yvars = 3; % state variables, one for obs variables
alpha = .05; %start with .05 maybe move down to .01
r1x = chi2inv(alpha/2,Nruns*xvars)/Nruns;
r2x = chi2inv(1-alpha/2,Nruns*xvars)/Nruns;
r1y = chi2inv(alpha/2, Nruns*yvars )./ Nruns;
r2y = chi2inv(1-alpha/2, Nruns*yvars )./ Nruns;
% Tuning:
QnewE = Qtrue*0.925;
RnewE = Rtrue*0.995;
QnewL = Qtrue*1;
RnewL = Rtrue*1.5;
P0E = diag([3 0.00001 1 0.00001])*1e-1;
P0L = diag([3 0.0001 1 0.0001])*1e-1;

%% NEES NIS Test for LKF and EKF
n = length(xstar);
Q = chol(Qtrue);
R = chol(Rtrue);
rng(100)
clear ytruth

for i = 1:Nruns
    % Step 1: Simulate Process Noise
    i
    tic
    xtruth(:,1) = xNLnoise(:,1);
    for j = 1:(n-1)
        noise = Q*randn(2,1); % 2x1 vector of noise, to be added to accel. terms
        tvec = [(j-1)*dt j*dt];
        %opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
        [~,xtruth_k] = ode45(@(t,x) OrbitNLnoise(t,x,mu,noise),tvec,xtruth(:,j));%,opts);
        xtruth(:,j+1) = xtruth_k(end,:)';
    % Step 2: Observations based on that new truth data
        t = (j-1)*dt;
        [~,~,~,ytruth{j+1}] = assembleH(t,xtruth(:,j+1));
    % Step 3: Add observation noise
        yNLnoise{j+1} = ytruth{j+1};
        % Add Observation Noise
        if numel(ytruth{j+1}) == 4
            yNLnoise{j+1}(1:3) = ytruth{j+1}(1:3) + R*randn(3,1);
        end
        if numel(ytruth{j+1}) > 4
            yNLnoise{j+1}(1:3,:) = ytruth{j+1}(1:3,:) + R*randn(3,2);
        end
    end
    init = xtruth(:,1)-xstar(:,1);
    [dx,sigmaLKF,~,NISL(i,:),PL] = KF(yNLnoise,ystar,init,F,P0L,H,length(xstar),QnewL,RnewL,Gamma,xtruth,xstar,ytruth,1);
    
    [xEKF,sigmaEKF,~,NISE(i,:),PE] = EKF(xtruth(:,1)+[1 0.1 1 0.1]'*.001,yNLnoise,P0E,QnewE,RnewE,Gamma,dt,mu,xtruth,1);
    % Step 5: Solve for NEES stuff
    exkL = xtruth - (xstar + dx); % NO clue about this one...
    for j = 1:(length(exkL))
        NEESL(i,j) = exkL(:,j)'*inv(PL{j})*exkL(:,j);
    end
    exkE = xtruth - xEKF;
    for j = 1:(length(exkE))
        NEESE(i,j) = exkE(:,j)'*inv(PE{j})*exkE(:,j);
    end
    toc
    % Repeat (mcs) number of times
end

%% Plot Results NEES NIS LKF
averageNEESL = mean(NEESL,1);
averageNISL = mean(NISL,1);
meanNEESL = mean(averageNEESL(1:50));
meanNISL = mean(averageNISL(1:50));
figure(3)
hold on
plot(1:length(xstar),averageNEESL,'*')
plot(1:length(xstar),ones(1,length(xstar))*r1x,'.')
plot(1:length(xstar),ones(1,length(xstar))*r2x,'.')
plot(1:length(xstar),ones(1,length(xstar))*meanNEESL,'-.k','LineWidth',1.5)
title('NEES Linearized Kalman Filter')

figure(4)
hold on
plot(2:length(xstar),averageNISL,'*')
plot(1:length(xstar),ones(1,length(xstar))*r1y,'.')
plot(1:length(xstar),ones(1,length(xstar))*r2y,'.')
plot(1:length(xstar),ones(1,length(xstar),1)*meanNISL,'-.k','LineWidth',1.5)
title('NIS Linearized Kalman Filter')

%% Plot Results NEES NIS EKF
averageNEESE = mean(NEESE,1);
averageNISE = mean(NISE,1);
meanNEESE = mean(averageNEESE);
meanNISE = mean(averageNISE);
figure(5)
hold on
plot(1:length(xstar),averageNEESE,'*')
plot(1:length(xstar),ones(1,length(xstar))*r1x,'.')
plot(1:length(xstar),ones(1,length(xstar))*r2x,'.')
plot(1:length(xstar),ones(1,length(xstar))*meanNEESE,'-.k','LineWidth',1.5)
title('NEES EKF')

figure(6)
hold on
plot(2:length(xstar),averageNISE,'*')
plot(1:length(xstar),ones(1,length(xstar))*r1y,'.')
plot(1:length(xstar),ones(1,length(xstar))*r2y,'.')
plot(1:length(xstar),ones(1,length(xstar),1)*meanNISE,'-.k','LineWidth',1.5)
title('NIS EKF')



%% Filter Runs:

%% Noisey Linearized Kalman Filter
% Function Call
[dxLKF, sigmaLKF, ~, ~] = KF(yNLnoise,ystar,dx(:,1),F,P0L,H,n,QnewL,RnewL,Gamma,xNL,xstar,ytruth,0);
xLKF = xstar + dxLKF;
% Plots:
figure(7) % Plot of kf output to Perturbed NL simulation error
subplot(4,1,1)
plot(tstar(2:end),1000*(xLKF(1,2:end)-xtruth(1,2:end)))
hold on
plot(tstar(2:end),0*(1:(length(tstar(2:end)))),'-.k')
hold on
plot(tstar(2:end),2000*sigmaLKF(1,2:end),'r')
plot(tstar(2:end),2000*(-sigmaLKF(1,2:end)),'r')
legend('X Distance LKF Error','2 Sigma Bounds')
title('Linearized Kalman Filter X Position Error')
xlabel('Time, [s]')
ylabel('Distance Error, [m]')
%ylim([-.25 .25])
subplot(4,1,2)
plot(tstar(2:end),1000*(xLKF(2,2:end)-xtruth(2,2:end)))
hold on
plot(tstar(2:end),2000*(sigmaLKF(2,2:end)),'r')
plot(tstar(2:end),2000*(-sigmaLKF(2,2:end)),'r')
legend('X Velocity LKF Error','2 Sigma Bounds')
title('Linearized Kalman Filter X Velocity Error')
xlabel('Time, [s]')
ylabel('Velocity Error, [m/s]')
%ylim([-.002 .002])
subplot(4,1,3)
plot(tstar(2:end),1000*(xLKF(3,2:end)-xtruth(3,2:end)))
hold on
plot(tstar(2:end),0*(1:(length(tstar(2:end)))),'-.k')
hold on
plot(tstar(2:end),2000*sigmaLKF(3,2:end),'r')
plot(tstar(2:end),2000*(-sigmaLKF(3,2:end)),'r')
legend('Y Distance LKF Error','2 Sigma Bounds')
title('Linearized Kalman Filter Y Position Error')
xlabel('Time, [s]')
ylabel('Distance Error, [m]')
%ylim([-.25 .25])
subplot(4,1,4)
plot(tstar(2:end),1000*(xLKF(4,2:end)-xtruth(4,2:end)))
hold on
plot(tstar(2:end),2000*(sigmaLKF(4,2:end)),'r')
plot(tstar(2:end),2000*(-sigmaLKF(4,2:end)),'r')
legend('Y Velocity LKF Error','2 Sigma Bounds')
title('Linearized Kalman Filter Y Velocity Error')
xlabel('Time, [s]')
ylabel('Velocity Error, [m/s]')
%ylim([-.002 .002])

%% Noisey EKF "Typical Situation"
% EKF Call
[xEKF,sigmaEKF,~,~,~] = EKF(xtruth(:,1)+[1 0.1 1 0.1]'*.001,yNLnoise,P0E,QnewE,RnewE,Gamma,dt,mu,xtruth,0);
% Plot Results
figure(8)
subplot(4,1,1)
plot(tstar(2:end),1000*(xEKF(1,2:end)-xtruth(1,2:end)))
hold on
plot(tstar(2:end),0*(1:(length(tstar(2:end)))),'-.k')
hold on
plot(tstar(2:end),2000*sigmaEKF(1,2:end),'r')
plot(tstar(2:end),2000*(-sigmaEKF(1,2:end)),'r')
legend('X Distance EKF Error','2 Sigma Bounds')
title('Extended Kalman Filter X Position Error')
xlabel('Time, [s]')
ylabel('Distance Error, [m]')
%ylim([-.25 .25])
subplot(4,1,2)
plot(tstar(2:end),1000*(xEKF(2,2:end)-xtruth(2,2:end)))
hold on
plot(tstar(2:end),2000*(sigmaEKF(2,2:end)),'r')
plot(tstar(2:end),2000*(-sigmaEKF(2,2:end)),'r')
legend('X Velocity EKF Error','2 Sigma Bounds')
title('Extended Kalman Filter X Velocity Error')
xlabel('Time, [s]')
ylabel('Velocity Error, [m/s]')
%ylim([-.002 .002])
subplot(4,1,3)
plot(tstar(2:end),1000*(xEKF(3,2:end)-xtruth(3,2:end)))
hold on
plot(tstar(2:end),0*(1:(length(tstar(2:end)))),'-.k')
hold on
plot(tstar(2:end),2000*sigmaEKF(3,2:end),'r')
plot(tstar(2:end),2000*(-sigmaEKF(3,2:end)),'r')
legend('Y Distance EKF Error','2 Sigma Bounds')
title('Extended Kalman Filter Y Position Error')
xlabel('Time, [s]')
ylabel('Distance Error, [m]')
%ylim([-.25 .25])
subplot(4,1,4)
plot(tstar(2:end),1000*(xEKF(4,2:end)-xtruth(4,2:end)))
hold on
plot(tstar(2:end),2000*(sigmaEKF(4,2:end)),'r')
plot(tstar(2:end),2000*(-sigmaEKF(4,2:end)),'r')
legend('Y Velocity EKF Error','2 Sigma Bounds')
title('Extended Kalman Filter Y Velocity Error')
xlabel('Time, [s]')
ylabel('Velocity Error, [m/s]')
%ylim([-.002 .002])


%% Provided Data: 

%% LKF
xvec = [Xnom dXnom Ynom dYnom]'; % X dX Y dY
% Call
t = 0;
n = length(ydata);
tvec = t:dt:dt*(length(ydata)-1); % time vector for nominal (NonLinear) simulation
opts = odeset('MaxStep',1e-2,'AbsTol',1e-12,'RelTol',1e-12);
[tstar,xstar] = ode45(@(t,x) OrbitNL(t,x,mu),tvec,xvec,opts);
xstar = xstar'; tstar = tstar'; % correct orientation, now many column vectors

% Create F Matrices: 
for k = 1:(n-1) % LINEARIZED DYNAMICS
    % Solve for A
    X = xstar(1,k);
    Y = xstar(3,k);
    A = [0 1 0 0;
        mu*(2*X^2-Y^2)/((Y^2+X^2).^(5/2)) 0 3*mu*X*Y/((X^2+Y^2).^(5/2)) 0;
        0 0 0 1;
        3*mu*X*Y/((X^2+Y^2).^(5/2)) 0 mu*(2*Y^2-X^2)/((Y^2+X^2).^(5/2)) 0];
    F{k} = eye(4) + dt*A;
    % Measurements
    t = (k-1)*10;
    [H{k+1},included{k+1},ystar_{k+1},ystar{k+1}] = assembleH(t,xstar(:,k+1));
end
% Call LKF
[dxLKF, sigmaLKF, ~, ~] = KF(ydata,ystar,dx(:,1)*15,F,P0L,H,n,QnewL,RnewL,Gamma,xNL,xstar,ytruth,0);
clear xLKF
xLKF(:,1) = xstar(:,1);
xLKF = xstar + dxLKF;

% Plot Results
figure(9)
subplot(4,1,1)
plot(dt*(1:length(xLKF)),xLKF(1,:))
hold on
title('Provided Data X Distance vs Time, LKF')
xlabel('Time, [s]')
ylabel('Distance, [km]')
subplot(4,1,2)
plot(dt*(1:length(xLKF)),xLKF(2,:))
hold on
title('Provided Data X Velocity vs Time, LKF')
xlabel('Time, [s]')
ylabel('Velocity, [km/s]')
subplot(4,1,3)
plot(dt*(1:length(xLKF)),xLKF(1,:))
hold on
title('Provided Data Y Distance vs Time, LKF')
xlabel('Time, [s]')
ylabel('Distance, [km]')
subplot(4,1,4)
plot(dt*(1:length(xLKF)),xLKF(2,:))
title('Provided Data Y Velocity vs Time, LKF')
xlabel('Time, [s]')
ylabel('Velocity, [km/s]')

figure(10)
subplot(4,1,1)
plot(dt*(1:length(xLKF)),2000*sigmaLKF(1,:))
hold on
title('Provided Data X Distance 2 Sigma Bound, LKF')
xlabel('Time, [s]')
ylabel('Distance, [m]')
subplot(4,1,2)
plot(dt*(1:length(xLKF)),2000*sigmaLKF(2,:))
hold on
title('Provided Data X Velocity 2 Sigma Bound, LKF')
xlabel('Time, [s]')
ylabel('Velocity, [m/s]')
subplot(4,1,3)
plot(dt*(1:length(xLKF)),2000*sigmaLKF(1,:))
hold on
title('Provided Data Y Distance 2 Sigma Bound, LKF')
xlabel('Time, [s]')
ylabel('Distance, [m]')
subplot(4,1,4)
plot(dt*(1:length(xLKF)),2000*sigmaLKF(2,:))
title('Provided Data Y Velocity 2 Sigma Bound, LKF')
xlabel('Time, [s]')
ylabel('Velocity, [m/s]')


%% EKF Call
[xEKF,sigmaEKF,~,~,~] = EKF(xNL(:,1),ydata,P0E,QnewE,RnewE,Gamma,dt,mu,xNL,0);
% Plot Results
figure(11)
subplot(4,1,1)
plot(dt*(1:length(xEKF)),xEKF(1,:))
hold on
title('Provided Data X Distance vs Time')
xlabel('Time, [s]')
ylabel('Distance, [km]')
subplot(4,1,2)
plot(dt*(1:length(xEKF)),xEKF(2,:))
hold on
title('Provided Data X Velocity vs Time')
xlabel('Time, [s]')
ylabel('Velocity, [km/s]')
subplot(4,1,3)
plot(dt*(1:length(xEKF)),xEKF(1,:))
hold on
title('Provided Data Y Distance vs Time')
xlabel('Time, [s]')
ylabel('Distance, [km]')
subplot(4,1,4)
plot(dt*(1:length(xEKF)),xEKF(2,:))
title('Provided Data Y Velocity vs Time')
xlabel('Time, [s]')
ylabel('Velocity, [km/s]')

figure(12)
subplot(4,1,1)
plot(dt*(1:length(xEKF)),2000*sigmaEKF(1,:))
hold on
title('Provided Data X Distance 2 Sigma Bound')
xlabel('Time, [s]')
ylabel('Distance, [m]')
subplot(4,1,2)
plot(dt*(1:length(xEKF)),2000*sigmaEKF(2,:))
hold on
title('Provided Data X Velocity 2 Sigma Bound')
xlabel('Time, [s]')
ylabel('Velocity, [m/s]')
subplot(4,1,3)
plot(dt*(1:length(xEKF)),2000*sigmaEKF(1,:))
hold on
title('Provided Data Y Distance 2 Sigma Bound')
xlabel('Time, [s]')
ylabel('Distance, [m]')
subplot(4,1,4)
plot(dt*(1:length(xEKF)),2000*sigmaEKF(2,:))
title('Provided Data Y Velocity 2 Sigma Bound')
xlabel('Time, [s]')
ylabel('Velocity, [m/s]')

