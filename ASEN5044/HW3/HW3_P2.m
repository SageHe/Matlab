%% Problem 2
clear all;close all;clc
J = 10;
F = 100;
T = 10;
x = [0;0];

A = [0 1;0 -F/J];
B = [0;1/J];

dx = 0.5;
for i = 1:(5/dx)
    xdot = A*x(:,i) + B*T;
    x(:,i+1) = x(:,i) + xdot*dx;
end
t = 0:dx:5;
figure
hold on
plot(t,x(1,:))
title('Angular Position (rads')
xlabel('Time (s)')
ylabel('Radians/ (Rad/S)')
plot(t,x(2,:))
title('Angular Rate VS Time \Deltat=0.5')
legend('Angular Pos.','Angular Rate')
%% Solving system using ss and lsim
% J = 10;
% F = 100;
% A = [0 1;0 -F/J];
% B = [0;1/J];
% C = [1 0;0 1];
% D = [0;0];
% T = 10;
% sys = ss(A,B,C,D);
% 
% dt = 0.005;
% 
% t = [0:dt:5];
% 
% u = ones(length(t),1)*T;
% 
% lsim(sys,u,t)

[V,D] = eig(A);

