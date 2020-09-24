%% CC 1, Q4
clear;close all;clc
addpath('C:\Users\sageh\OneDrive\Documents\MATLAB\ASEN5010\Matlab')
sigma_BN0 = [0.1 0.2 -0.1]';
w_BN0 = deg2rad([30 10 -20]');
gains = [10 5]'; %gains in order [P K]'
tstart = 0;
tend = 120;
[X,~] = RK4_HW6(sigma_BN0,w_BN0,tstart,tend,.1,gains);

figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))

close all;clc

[X,~,err] = RK4_HW6(sigma_BN0,w_BN0,tstart,tend,.01,gains);
norm(err(:,1,4001));

close all;clc

[X,~,err] = RK4_HW6(sigma_BN0,w_BN0,tstart,tend,.01,gains);
norm(err(:,1,2001));

clc

[X,~,err] = RK4_HW6(sigma_BN0,w_BN0,tstart,tend,.01,gains);
norm(err(:,1,8001));

clc
[X,~,err] = RK4_HW6(sigma_BN0,w_BN0,tstart,tend,.01,gains);
norm(err(:,1,3001))