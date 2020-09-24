%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 4: BALANCED AND UNBALANCED WHEEL
%
% Created:  02/28/2018 - Charles Puskar
% Modified: 03/05/2018 - Charles Puskar
%
% lab4.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% PLOTTING MODEL
theta = linspace(0,15,1000);

figure(1)
hold on
grid minor
xlabel('Angular Position [rad]')
ylabel('Expected Angular Velocity [rad/s]')
plot(theta,model_1(theta),theta,model_2(theta),theta,model_3(theta),theta,model_4(theta),'LineWidth',1);
legend({'model\_1','model\_2','model\_3','model\_4'},'Location','Southeast')


%% READING DATA
file_bal1 = 'balanced1';
file_bal2 = 'balanced2';
file_unbal1 = 'unbalanced1';
file_unbal2 = 'unbalanced2';

[theta_bal1, w_bal1] = read_data(file_bal1);
[theta_bal2, w_bal2] = read_data(file_bal2);
[theta_unbal1, w_unbal1] = read_data(file_unbal1);
[theta_unbal2, w_unbal2] = read_data(file_unbal2);

%% PLOTTING BALANCED DATA
figure(2)
hold on
grid minor
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
plot(theta_bal1,w_bal1,'*',theta_bal2,w_bal2,'*','LineWidth',1)
plot(theta,model_1(theta),'--',theta,model_2(theta),'-','LineWidth',1)
legend({'trial 1 data','trial 2 data','model\_1','model\_2'},'Location','Southeast')

%% PLOTTING UNBALANCED DATA
figure(3)
hold on
grid minor
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
plot(theta_unbal1,w_unbal1,'*',theta_unbal2,w_unbal2,'*','LineWidth',1)
plot(theta,model_3(theta),'--',theta,model_4(theta),'-','LineWidth',1)
legend({'trial 1 data','trial 2 data','model\_3','model\_4'},'Location','Southeast')

%% RESIDUALS

