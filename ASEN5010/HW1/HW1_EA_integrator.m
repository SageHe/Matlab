clear all;close all;clc

Eulers = [40 30 80]';
Eulers = deg2rad(Eulers);
N = 10000;
t = linspace(0,42,N);
deltat = t(2) - t(1);
for i = 1:N
    Edots = 1/cos(Eulers(2,i))*[0 sin(Eulers(3,i)) cos(Eulers(3,i));...
                            0 cos(Eulers(3,i))*cos(Eulers(2,i)) -sin(Eulers(3,i))*cos(Eulers(2,i));...
                            cos(Eulers(2,i)) sin(Eulers(3,i))*sin(Eulers(2,i)) cos(Eulers(3,i))*sin(Eulers(2,i))];
    omega = [sin(0.1*t(i)) 0.01 cos(0.1*t(i))]*deg2rad(20);
    omega = omega';
    
    Edots = Edots*omega;
    
    Eulers(:,(i+1)) = Eulers(:,i) + Edots*deltat;
end
% plot(t,Eulers)

x = Eulers(:,end);

% x = deg2rad(x);

norm(x)