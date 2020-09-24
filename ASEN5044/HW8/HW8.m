%% HW8 part a
clear;close all;clc

Aa = [0 1 0 0;0 0 0 -0.045;0 0 0 1;0 0.045 0 0];
Ab = [0 1 0 0;0 0 0 0.045;0 0 0 1;0 -0.045 0 0];

gamma = [0 0;1 0;0 0;0 1];

W = 10*[2 0.05;0.05 0.5];

dt = 0.5;

Za = [-Aa gamma*W*gamma';zeros(4) Aa'];
Za = dt*Za;

Zb = [-Ab gamma*W*gamma';zeros(4) Ab'];
Zb = dt*Zb;

eza = expm(Za);
ezb = expm(Zb);

Fa = eza([5:8],[5:8])';
Qa = eza([5:8],[5:8])'*eza([1:4],[5:8]);
Fb = ezb([5:8],[5:8])';
Qb = ezb([5:8],[5:8])'*ezb([1:4],[5:8]);
%% Part bi
%load in data
data = load('hw8Problem1_data.mat');
% fix random number seed to 100
rng(100)
R = [20 0.05;0.05 20];
Sv = chol(R,'lower');
H = [1 0 0 0;0 0 1 0];
for i = 1:201
    y(:,i) = H*data.xasingle_truth(:,i) + Sv*randn(2,1);
end
t = [0:.5:20];
figure
plot(t,y(:,[1:41]))
xlabel('Time (s)')
ylabel('Position (m)')
title('Position Measurement VS Time')
legend('North Mesurement','East Measurement')
%% Part bii
mu0 = [0 85*cos(pi/4) 0 -85*sin(pi/4)]';
P0 = 900*diag([10 2 10 2]);
xhat = mu0;
Pkm1p = P0;
sigmas = [P0(1,1) P0(2,2) P0(3,3) P0(4,4)];
for i = 2:201
    Pkm = Fa*Pkm1p*Fa' + Qa;
    Kk = Pkm*H'*inv(H*Pkm*H'+R);
    Xkm = Fa*xhat(:,i-1);
    Xkp = Xkm + Kk*(y(:,i) - H*Xkm);
    Pkp = (eye(4) - Kk*H)*Pkm;
    sigmas = [sigmas;Pkp(1,1) Pkp(2,2) Pkp(3,3) Pkp(4,4)];
    xhat(:,i) = Xkp;
    Pkm1p = Pkp;
end
ekp = data.xasingle_truth - xhat;
t = [0:0.5:100];
figure
subplot(4,1,1)
hold on
plot(t,ekp(1,:))
plot(t,2*sqrt(sigmas(:,1)),'--')
plot(t,-2*sqrt(sigmas(:,1)),'--')
xlabel('Time (s)')
ylabel('East Pos. (m)')
title('Position and Velocity Error VS Time')
ylim([-25 25])
subplot(4,1,2)
hold on
plot(t,ekp(2,:))
plot(t,2*sqrt(sigmas(:,2)),'--')
plot(t,-2*sqrt(sigmas(:,2)),'--')
xlabel('Time (s)')
ylabel('East Vel. (m/s)')
ylim([-25 25])
subplot(4,1,3)
hold on
plot(t,ekp(3,:))
plot(t,2*sqrt(sigmas(:,3)),'--')
plot(t,-2*sqrt(sigmas(:,3)),'--')
xlabel('Time (s)')
ylabel('North Pos. (m)')
ylim([-25 25])
subplot(4,1,4)
hold on
plot(t,ekp(4,:))
plot(t,2*sqrt(sigmas(:,4)),'--')
plot(t,-2*sqrt(sigmas(:,4)),'--')
xlabel('Time (s)')
ylabel('North Vel. (m/s)')
ylim([-25 25])
%% Part ci
% define givens
muA0 = [0 85*cos(pi/4) 0 -85*sin(pi/4)]';
PA0 = 900*diag([10 2 10 2]);
muB0 = [3200 85*cos(pi/4) 3200 -85*sin(pi/4)]';
PB0 = 900*diag([11 4 11 4]);
Rd = [10 0.15;0.15 10];
Sv = chol(Rd,'lower');
H = [1 0 0 0;0 0 1 0];
for i = 1:201
    yd(:,i) = H*data.xadouble_truth(:,i) - H*data.xbdouble_truth(:,i) + Sv*randn(2,1);
    yAp(:,i) = H*data.xadouble_truth(:,i) + Sv*randn(2,1);
end
ys = [yAp;yd];
% Augment new F,Q,H,R, and P(0) matrices based off of existing matries 
F = blkdiag(Fa,Fb);
Q = blkdiag(Qa,Qb);
H = [1 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0;1 0 0 0 -1 0 0 0;0 0 1 0 0 0 -1 0];
R = blkdiag(R,Rd);
P0 = blkdiag(PA0,PB0);
mu0 = [muA0;muB0];
xshat = mu0;
Pkm1p = P0;
for i = 2:201
    Pkm = F*Pkm1p*F' + Q;
    Kk = Pkm*H'*inv(H*Pkm*H'+R);
    Xkm = F*xshat(:,i-1);
    Xkp = Xkm + Kk*(ys(:,i) - H*Xkm);
    Pkp = (eye(8) - Kk*H)*Pkm;
    xshat(:,i) = Xkp;
    Pkm1p = Pkp;
end
ekp_A = [data.xadouble_truth - xshat([1:4],:)];
ekp_B = [data.xbdouble_truth - xshat([5:8],:)];
t = [0:0.5:100];
figure
hold on
plot(t,ekp_A(1,:))
plot(t,ekp_A(3,:))
ylim([-25 25])
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Position Error VS Time (s), Plane A')
legend('East Pos. Error','North Pos. Error')
figure
hold on
plot(t,ekp_B(1,:))
plot(t,ekp_B(3,:))
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Position Error VS Time (s), Plane B')
legend('East pos, Error','North Pos. Error')
ylim([-25 25])
%% part cii
muA0 = [0 85*cos(pi/4) 0 -85*sin(pi/4)]';
PA0 = 900*diag([10 2 10 2]);
muB0 = [3200 85*cos(pi/4) 3200 -85*sin(pi/4)]';
PB0 = 900*diag([11 4 11 4]);
Rd = [10 0.15;0.15 10];
Sv = chol(Rd,'lower');
H = [1 0 0 0;0 0 1 0];
for i = 1:201
    yd(:,i) = H*data.xadouble_truth(:,i) - H*data.xbdouble_truth(:,i) + Sv*randn(2,1);
end
F = blkdiag(Fa,Fb);
Q = blkdiag(Qa,Qb);
H = [1 0 0 0 -1 0 0 0;0 0 1 0 0 0 -1 0];
R = Rd;
P0 = blkdiag(PA0,PB0);
mu0 = [muA0;muB0];
xshat = mu0;
Pkm1p = P0;
for i = 2:201
    Pkm = F*Pkm1p*F' + Q;
    Kk = Pkm*H'*inv(H*Pkm*H'+R);
    Xkm = F*xshat(:,i-1);
    Xkp = Xkm + Kk*(yd(:,i) - H*Xkm);
    Pkp = (eye(8) - Kk*H)*Pkm;
    xshat(:,i) = Xkp;
    Pkm1p = Pkp;
end
ekp_A = [data.xadouble_truth - xshat([1:4],:)];
ekp_B = [data.xbdouble_truth - xshat([5:8],:)];
t = [0:0.5:100];
figure
hold on
plot(t,ekp_A(1,:))
plot(t,ekp_A(3,:))
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Position Error (m), Transponder Only, Plane A')
legend('East Pos. Error','North Pos. Error')
% ylim([-25 25])
figure
hold on
plot(t,ekp_B(1,:))
plot(t,ekp_B(3,:))
xlabel('Time (s)')
ylabel('Position Error (m)')
title('Position Error (m), Transponder Only, Plane B')
legend('East Pos. Error','North Pos. Error')
% ylim([-25 25])