clear all; close all;clc

Ix = 500; %Kg m^3
Iy = 750; %Kg m^3
Iz = 1000;%Kg m^3
p0 = 20; %rad/s
dt = 0.1; %s

A = [0 0 0;0 0 (p0*(Ix - Iz)/Iy); 0 (p0*(Iy - Ix)/Iz) 0]; %Calculate A matrix based on givene equations 

x0 = [0;0.1;0]; %Specify ICs
count = 1;
for i = 0:dt:5 %Continuously update STM and calculate x based on current t
    STM = expm(A*i);
    x(count,:) = STM*x0;
    count = count + 1;
end

t = [0:0.1:5];
%Plot results 
% figure
% subplot(1,3,1)
% plot(t,x(:,1))
% title('\Delta p VS Time')
% xlabel('Time (s)')
% ylabel('Rotation Rate (rad/s)')
% subplot(1,3,2)
% plot(t,x(:,2))
% title('\Delta q VS Time')
% xlabel('Time (s)')
% ylabel('Rotation Rate (rad/s)')
% subplot(1,3,3)
% plot(t,x(:,3))
% title('\Delta r VS Time')
% xlabel('Time (s)')
% ylabel('Rotation Rate (rad/s)')
figure
hold on
plot(t,x(:,1))
plot(t,x(:,2))
plot(t,x(:,3))
title('Rotation Rates Vs Time')
xlabel('Time (s)')
ylabel('Rotation Rate (rad/s)')
legend('\Delta p','\Delta q','\Delta r')