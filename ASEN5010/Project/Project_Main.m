%ASEN 5010 PROJECT - Main Script
% Author: Sage Herrin
%% Task 1
clear;close all;clc
addpath('C:\Users\sageh\OneDrive\Documents\MATLAB\ASEN5010\Matlab')
%Derive inertial spacecraft velocity vector rdot. Note for circular orbit
%thetadot is constant
global rLMO rGMO omegaLMO iLMO thetaLMO omegaGMO iGMO thetaGMO thetadotLMO thetadotGMO
rLMO = 3396.19 + 400;
omegaLMO = deg2rad(20);
iLMO = deg2rad(30);
thetaLMO = deg2rad(60);
thetadotLMO = 0.000884797;
thetadotGMO = 0.0000709003;
t = 450;
[rLM,rdotLM] = velandpos(rLMO,omegaLMO,iLMO,thetaLMO+thetadotLMO*t);
fid = fopen('N_Pos_450.txt','w');
fprintf(fid,'%f ',rLM');
fclose(fid);
fid = fopen('N_Vel_450.txt','w');
fprintf(fid,'%f ',rdotLM');
fclose(fid);
rGMO = 3396.19 + 17028.01;
thetadotGMO = 0.0000709003; %rad/s
omegaGMO = 0;
iGMO = 0;
thetaGMO = deg2rad(250);
t = 1150;
[rGM,rdotGM] = velandpos(rGMO,omegaGMO,iGMO,thetaGMO+thetadotGMO*t);
fid = fopen('N_Pos_1150.txt','w');
fprintf(fid,'%f ',rGM');
fclose(fid);
fid = fopen('N_Vel_1150.txt','w');
fprintf(fid,'%f ',rdotGM');
fclose(fid);
%% Task 2
% Determin the DCM HN as only a function of time 
HN = calcHN(300);
fid = fopen('HN.txt','w');
fprintf(fid,'%f ',HN');
fclose(fid);
%% Task 3
% Determine sun pointing reference frame Rs by defining the DCM [RsN]
%Given definition of R frame relative to N frame, DCM RsN can be easily
%determined as follows
RsN = [-1 0 0;0 0 1;0 1 0];
%Since Rs reference frame is fixed wrt to N frame, N frame derivitive of
%wRs/N is 0
wRsN = [0 0 0];
fid = fopen('RsN.txt','w');
fprintf(fid,'%f ',RsN');
fclose(fid);
fid = fopen('wRsN.txt','w');
fprintf(fid,'%f ',wRsN);
fclose(fid);
%% Task 4
%Determine nadir pointing reference frame Rn by defining DCM [RnN]
%A simple RH matrix can be deterined based on the description of the
%requirements, yielding
RH = [-1 0 0;0 1 0;0 0 -1];
%RnN can be determined via RnN = RH*HN, where HN can be determined from
%task 2
%resest thetaLMO for this scenario
% thetaLMO = deg2rad(60);
RnN = calcRnN(330);
fid = fopen('RnN.txt','w');
fprintf(fid,'%f ',RnN');
fclose(fid);
% thetaLMO = deg2rad(60);
NwRnN = calcNwRnN(330);
fid = fopen('NwRnN.txt','w');
fprintf(fid,'%.12f ',NwRnN');
fclose(fid);
%% Task 5
%Determine communication mode reference frame Rc by defining DCM [RcN]
t = 330;
RcN = calcRcN(t);
fid = fopen('RcN.txt','w');
fprintf(fid,'%.12f ',RcN');
fclose(fid);
NwRcN = calcNwRcN(t);
% determine angular velocity vector NwRcN
% RcN331 = calcRcN(331);
% RcN330 = calcRcN(330);
% 
% RcNdot = RcN331-RcN330;
% omegatilde = RcNdot*RcN330';
% NwRcN = [omegatilde(2,3) -omegatilde(1,3) omegatilde(1,2)]';
fid = fopen('NwRcN.txt','w');
fprintf(fid,'%.12f ',NwRcN');
fclose(fid);
%% Task 6
% compute attitude error by computing sigma and w for several frames
t = 0;
sigma_BN = [0.3 -0.4 0.5]';
Bw_BN = deg2rad([1.00 1.75 -2.20]');
[sigma_BR_s,w_BR_s] = att_err(sigma_BN,Bw_BN,RsN,wRsN');
fid = fopen('sigma_BR_s.txt','w');
fprintf(fid,'%.12f ',sigma_BR_s');
fclose(fid);
fid = fopen('w_BR_s.txt','w');
fprintf(fid,'%.12f ',w_BR_s');
fclose(fid);
RnN = calcRnN(0);
NwRnN = calcNwRnN(0);
[sigma_BR_n,w_BR_n] = att_err(sigma_BN,Bw_BN,RnN,NwRnN);
sigma_BR_n = -sigma_BR_n/(sigma_BR_n'*sigma_BR_n);
fid = fopen('sigma_BR_n.txt','w');
fprintf(fid,'%.12f ',sigma_BR_n');
fclose(fid)
fid = fopen('w_BR_n.txt','w');
fprintf(fid,'%.12f ',w_BR_n');
fclose(fid)
RcN = calcRcN(0);
NwRcN = calcNwRcN(0);
[sigma_BR_c,w_BR_c] = att_err(sigma_BN,Bw_BN,RcN,NwRcN);
fid = fopen('sigma_BR_c.txt','w');
fprintf(fid,'%.12f ',sigma_BR_c');
fclose(fid)
fid = fopen('w_BR_c.txt','w');
fprintf(fid,'%.12f ',w_BR_c');
fclose(fid);
%% Task 7 - RK4 integrator 
B_I = [10 0 0;0 5 0;0 0 7.5];
sigma_BN0 = [0.3 -0.4 0.5]';
Bw_BN0 = deg2rad([1.00 1.75 -2.20]');
tstart = 0;
tend = 500;
dt = 1;
u = [0 0 0]';
[X,u] = RK4(sigma_BN0,Bw_BN0,tstart,tend,dt,u);
figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))
xlabel('Time (s)')
ylabel('\omega_{B/N}')
legend('\omega_1','\omega_2','\omega_3')
title('Spacecraft Angular Rate VS Time')
sigma_BN = X(:,1,end);
Bw_BN = X(:,2,end);
H = B_I*Bw_BN;
fid = fopen('B_H.txt','w');
fprintf(fid,'%.12f ',H');
fclose(fid);
T = .5*Bw_BN'*B_I*Bw_BN;
fid = fopen('T.txt','w');
fprintf(fid,'%.12f ',T');
fclose(fid);
fid = fopen('sigma_BN.txt','w');
fprintf(fid,'%.12f ',sigma_BN');
fclose(fid);
BN = MRP2C(sigma_BN);
N_H = BN'*H;
fid = fopen('N_H.txt','w');
fprintf(fid,'%.12f ',N_H');
fclose(fid);
u = [0.01 -0.01 0.02]';
tend = 100;
[X,~] = RK4(sigma_BN0,Bw_BN0,tstart,tend,dt,u);
sigma_BN_u = X(:,1,end);
fid = fopen('sigma_BN_u.txt','w');
fprintf(fid,'%.12f ',sigma_BN_u');
fclose(fid);
%% Task 8 - Sun Pointing Control 
% first compute gains K and P using given info; all decay time constants
% need to be 120 seconds or less, all closed loop response for all sigma_BN
% components should be either crit. damped or under damped
T = 120;
zeta = 1;
% for i = 1:3
%     K(i) = ((2*sqrt(B_I(i,i)))/(T*zeta))^2;
% end
% K = max(K);
for i = 1:3
    P(i) = (2*B_I(i,i))/T;
    for j = 1:3
        Ti(j) = 2*B_I(j,j)/P(i);
    end
    if numel(Ti <= 120) == sum(Ti <= 120)
        P = P(i);
        break
    end
end
for i = 1:3
    K(i) = (P/(zeta*sqrt(B_I(i,i))))^2;
end
K = max(K);
% P = P';
gains = [P;K];
fid = fopen('Part8_GAINS.txt','w');
fprintf(fid,'%.12f ',gains');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,15,1,gains);
sigma_BN_15 = X(:,1,end);
fid = fopen('Tast8_sBN_15.txt','w');
fprintf(fid,'%.12f ',sigma_BN_15');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,100,1,gains);
sigma_BN_100 = X(:,1,end);
fid = fopen('Tast8_sBN_100.txt','w');
fprintf(fid,'%.12f ',sigma_BN_100');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,200,1,gains);
sigma_BN_200 = X(:,1,end);
fid = fopen('Tast8_sBN_200.txt','w');
fprintf(fid,'%.12f ',sigma_BN_200');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,400,1,gains);
sigma_BN_400 = X(:,1,end);
fid = fopen('Tast8_sBN_400.txt','w');
fprintf(fid,'%.12f ',sigma_BN_400');
fclose(fid);
figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))
xlabel('Time (s)')
ylabel('\omega_{B/N}')
title('Spacecraft Angular Rate Vs Time')
legend('\omega_1','\omega_2','\omega_3')
%% Task 9 - Nadir pointing control mode
[X,~] = RK4(sigma_BN0,Bw_BN0,0,15,1,gains);
sigma_BN_15 = X(:,1,end);
fid = fopen('Task9_nBN_15.txt','w');
fprintf(fid,'%.12f ',sigma_BN_15');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,100,1,gains);
sigma_BN_100 = X(:,1,end);
fid = fopen('Task9_nBN_100.txt','w');
fprintf(fid,'%.12f ',sigma_BN_100');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,200,1,gains);
sigma_BN_200 = X(:,1,end);
fid = fopen('Task9_nBN_200.txt','w');
fprintf(fid,'%.12f ',sigma_BN_200');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,400,1,gains);
sigma_BN_400 = X(:,1,end);
fid = fopen('Task9_nBN_400.txt','w');
fprintf(fid,'%.12f ',sigma_BN_400');
fclose(fid);
figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))
xlabel('Time (s)')
ylabel('\omega_{B/N}')
title('Spacecraft Angular Rate Vs Time')
legend('\omega_1','\omega_2','\omega_3')
% Task 10 - GMO pointing control mode
[X,~] = RK4(sigma_BN0,Bw_BN0,0,15,1,gains);
sigma_BN_15 = X(:,1,end);
fid = fopen('Task10_cBN_15.txt','w');
fprintf(fid,'%.12f ',sigma_BN_15');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,100,1,gains);
sigma_BN_100 = X(:,1,end);
fid = fopen('Task10_cBN_100.txt','w');
fprintf(fid,'%.12f ',sigma_BN_100');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,200,1,gains);
sigma_BN_200 = X(:,1,end);
fid = fopen('Task10_cBN_200.txt','w');
fprintf(fid,'%.12f ',sigma_BN_200');
fclose(fid);
[X,~] = RK4(sigma_BN0,Bw_BN0,0,400,1,gains);
sigma_BN_400 = X(:,1,end);
fid = fopen('Task10_cBN_400.txt','w');
fprintf(fid,'%.12f ',sigma_BN_400');
fclose(fid);
figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))
xlabel('Time(s)')
ylabel('\omega_{B/N}')
title('Spacecraft Angular Rate Vs Time')
legend('\omega_1','\omega_2','\omega_3')
%% Task 11 - 
[X,~] = RK4(sigma_BN0,Bw_BN0,0,6500,1,gains);
fid = fopen('Task11_sigma_BN_300.txt','w');
fprintf(fid,'%.12f ',X(:,1,301)');
fclose(fid)
fid = fopen('Task11_sigma_BN_2100.txt','w');
fprintf(fid,'%.12f ',X(:,1,2101)');
fclose(fid)
fid = fopen('Task11_sigma_BN_3400.txt','w');
fprintf(fid,'%.12f ',X(:,1,3401)');
fclose(fid)
fid = fopen('Task11_sigma_BN_4400.txt','w');
fprintf(fid,'%.12f ',X(:,1,4401)');
fclose(fid)
fid = fopen('Task11_sigma_BN_5600.txt','w');
fprintf(fid,'%.12f ',X(:,1,5601)');
fclose(fid)
figure
plot(squeeze(X(1,2,:)))
hold on
plot(squeeze(X(2,2,:)))
plot(squeeze(X(3,2,:)))
xlabel('Time (s)')
ylabel('\omega_{B/N}')
title('Spacecraft Angular Rate Vs Time')
legend('\omega_1','\omega_2','\omega_3')