%% part b
clear;close all;clc
mu = 398600;
r = 6678;
A = [0 1 0 0;-mu/(r^3) 0 0 0;0 0 0 1;0 0 -mu/(r^3) 0];
B = [0 0;1 0;0 0;0 1];
L = B;

dt = 10; %secs

Ahat = [A B;zeros(2,6)];
eAhat = expm(Ahat*dt);
F = eAhat([1:4],[1:4]);
G = eAhat([1:4],[5:6]);

[V,D] = eig(F);

P = [G F*G F^2*G F^3*G];

rank(P)
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
X = [X0];
%     for j = 1:12
%         theta0(j) = (j - 1)*(pi/6);
%         Xs(j) = Re*cos(We + theta0(j));
%         Xdots(j) = -Re*sin(We+theta0(j))*We;
%         Ys(j) = Re*sin(We + theta0(j));
%         Ydots(j) = Re*cos(We + theta0(j)); 
%         rho(j) = sqrt((X0(1)-Xs0(j))^2 + X0(3)-Ys0(j));
%         rhodot(j) = ((X0(1) - Xs0(j))*(X0(2)-Xdots0(j)) + (X0(3)-Ys0(j))*(X0(4)-Ydots0(j)))/rho0(j);
%         phi(j) = atan2((X0(3)-Ys0(j)),(X0(1)-Xs0(j)));
%     %     H{j} = [
%     end
%Determine H matrix and calculate observations y over one orbit period
y = [];
for i = 1:581
    X(:,i+1) = F^i*X0;
    yloop = [];
    for j = 1:12
        theta0(j) = (j - 1)*(pi/6);
        Xs = Re*cos(We*(i*10) + theta0(j));
        Xdots = -Re*sin(We*(i*10)+theta0(j))*We;
        Ys = Re*sin(We*(i*10) + theta0(j));
        Ydots = Re*cos(We*(i*10) + theta0(j)); 
        rho = sqrt((X(1,i)-Xs)^2 + (X(3,i)-Ys)^2);
        rhodot = ((X(1,i) - Xs)*(X(2,i)-Xdots) + (X(3,i)-Ys)*(X(4,i)-Ydots))/rho;
        num = ((X(1,i) - Xs)*(X(2,i)-Xdots) + (X(3,i)-Ys)*(X(4,i)-Ydots));
        phi = atan2((X(3,i)-Ys),(X(1,i)-Xs));
        H = [(X(1,i)-Xs)/rho 0 (X(3,i)-Ys)/rho 0;
                (rho*(X(2,i)-Xdots)-(num)*((X(1,i)-Xs)/rho))/(rho^2) (X(1,i)-Xs)/rho  (rho*(X(4,i)-Ydots)-(num)*((X(3,i)-Ys)/rho))/(rho^2) (X(3,i)-Ys)/(rho^2);
                (Ys-X(3,i))/(rho^2) 0 (X(1,i)-Xs)/(rho^2) 0];
        if phi >= ((-pi/2)+atan2(Ys,Xs)) && phi <= ((pi/2)+atan2(Ys,Xs))
            ytemp = H*(X(:,i) - [Xs Xdots Ys Ydots]');
        else
            ytemp = nan(3,1);
        end
        yloop = [yloop; ytemp];
    end
    y = [y yloop];
end
tvec = linspace(0,581*10,581);
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