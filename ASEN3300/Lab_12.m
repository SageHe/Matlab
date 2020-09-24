clear all;close all;clc
%call gps_newparser function to interpret GPS text files
DLC = GPS_newparser('GPSLOG01.txt','out1.txt');
walk = GPS_newparser('Main_walk_ec_loop.txt','out2.txt');
%convert lat,long, and height into x,y,z
dlcxyz = lla2ecef([DLC.lat;DLC.long;DLC.alt]');
walkxyz = lla2ecef([walk.lat;walk.long;walk.alt]');
time = DLC.time';
%determine mean and standard dev. of x,y,and z
xmean = mean(dlcxyz(:,1));
xsigma = std(dlcxyz(:,1));
ymean = mean(dlcxyz(:,2));
ysigma = std(dlcxyz(:,2));
zmean = mean(dlcxyz(:,3));
zsigma = std(dlcxyz(:,3));
%pulling from reference positions after post processing
x_ref = ones(1,length(time))*-1288106.079;
y_ref = ones(1,length(time))*-4720833.680;
z_ref = ones(1,length(time))*4079649.515;
%plot static x,y,z vs time marble ball
figure(1)
hold on
plot(time,dlcxyz(:,1))
plot(time,x_ref)
title('Marbe Ball X, mean = -1.2881e6, \sigma = 4.009')
xlabel('Time (s)')
ylabel('X (meters)')
legend('Data','Reference')
figure(2)
hold on
plot(time,dlcxyz(:,2))
plot(time,y_ref)
title('Marble Ball Y,mean = -4.7209e6, \sigma = 2.0396')
xlabel('Time (s)')
ylabel('Y (meters)')
legend('Data','Reference')
figure(3)
hold on
plot(time,dlcxyz(:,3))
plot(time,z_ref)
title('Mable Ball Z, mean = 4.0797e6, \sigma = 1.8433')
xlabel('Time (s)')
ylabel('Z (meters)')
legend('Data','Legend')
%data for other group is GPSLOG08 file
BL = GPS_newparser('GPSLOG08.txt','out1.txt');
BLxyz = lla2ecef([BL.lat;BL.long;BL.alt]');
time = BL.time;
%determine mean and standard dev. of x,y,and z
xmean = mean(BLxyz(:,1));
xsigma = std(BLxyz(:,1));
ymean = mean(BLxyz(:,2));
ysigma = std(BLxyz(:,2));
zmean = mean(BLxyz(:,3));
zsigma = std(BLxyz(:,3));
%pulling from reference positions after post processing
x_ref = ones(1,length(time))*-1288198.742;
y_ref = ones(1,length(time))*-4721358.952;
z_ref = ones(1,length(time))*4079036.732;
%plot static x,y,z vs time marble ball
figure(4)
hold on
plot(time,BLxyz(:,1))
plot(time,x_ref)
% ylim([-1288200 -1288100])
title('Baseline X, mean = -1.2882e6, \sigma = 0.5059')
xlabel('Time (s)')
ylabel('X (meters)')
legend('Data','Reference')
figure(5)
hold on
plot(time,BLxyz(:,2))
plot(time,y_ref)
title('Baseline Y,mean = -4.7214e6, \sigma = 1.4955')
xlabel('Time (s)')
ylabel('Y (meters)')
legend('Data','Reference')
figure(6)
hold on
plot(time,BLxyz(:,3))
plot(time,z_ref)
title('Baseline Z, mean = 4.0791e6, \sigma = 0.3821')
xlabel('Time (s)')
ylabel('Z (meters)')
legend('Data','Reference')