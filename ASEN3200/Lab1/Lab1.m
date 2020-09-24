%Question 5.1, measuring bias and sensitivity of MEMS gyro
clear all; close all; clc
Run1 = load('AM_Unit06_MEMS_COMP_RUN1');
Run1([1:31],:) = [];
time = Run1(:,1);
gyro_output = -Run1(:,2); %rads/s
input_rate = Run1(:,3); %rpm

input_rate = input_rate*(1/60)*(2*pi); %convert rpm to rads/s
figure(1)
subplot(3,1,1)
hold on
plot(input_rate,gyro_output)
title('Input rate VS Gyro output, .5A .2Hz')
xlabel('Encoder Input rad/ms')
ylabel('Gyro Output rad/ms')

P = polyfit(input_rate,gyro_output,1);
y = polyval(P,input_rate);
bias = [P(2)];
sens = [P(1)];
plot(input_rate,y,'LineWidth',2)
legend('Data','Best fit line')
adjust1 = gyro_output*(1/sens) - bias;
figure(2)
subplot(3,1,1)
hold on
plot(Run1(:,1),input_rate)
plot(Run1(:,1),gyro_output)
plot(Run1(:,1),adjust1)
title('Encoder Vs Gyro Vs Adjusted Angular Rates, .5A .2Hz')
xlabel('Time (ms)')
ylabel('Angular Velocity (rad/ms)')
legend('Input Rate','Output Rate','Adjusted')
for i = 2:length(time)
    posG1(i) = trapz([time(1:i)],[gyro_output(1:i)]);
    posE1(i) = trapz([time(1:i)],[input_rate(1:i)]);
    adjpos1(i) = trapz([time(1:i)],[adjust1(1:i)]);
end
figure(3)
subplot(3,1,1)
hold on
plot(time,posG1)
plot(time,posE1)
plot(time,adjpos1)
title('Angular Position Vs Time, .5A .2Hz')
xlabel('Time (ms)')
ylabel('Radians')
legend('Gyro Output','Encoder Input','Adjusted Gyro Readings')
%Second run testing mems gyro
Run2 = load('AM_Unit06_MEMS_COMP_RUN02');
Run2(1,:) = [];
time = Run2(:,1);
gyro_output = -Run2(:,2); %rads/s
input_rate = Run2(:,3); %rpm

input_rate = input_rate*(1/60)*(2*pi); %convert rpm to rads/s
figure(1)
subplot(3,1,2)
hold on
plot(input_rate,gyro_output)
title('Input rate VS Gyro output, .75A .2Hz')
xlabel('Encoder Input rad/ms')
ylabel('Gyro Output rad/ms')
legend('Data','Fit line')

P = polyfit(input_rate,gyro_output,1);
y = polyval(P,input_rate);
bias = [bias P(2)];
sens = [sens P(1)];
plot(input_rate,y,'LineWidth',2)
legend('Data','Best fit line')

adjust2 = gyro_output*(1/sens(end)) - bias(end);
figure(2)
subplot(3,1,2)
hold on
plot(Run2(:,1),input_rate)
plot(Run2(:,1),gyro_output)
plot(Run2(:,1),adjust2)
title('Encoder Vs Gyro Vs Adjusted Angular Rates, .7A .2Hz')
xlabel('Time (ms)')
ylabel('Angular Velocity (rad/ms)')
legend('Input Rate','Gyro Output','Adjusted')
for i = 2:length(time)
    posG2(i) = trapz([time(1:i)],[gyro_output(1:i)]);
    posE2(i) = trapz([time(1:i)],[input_rate(1:i)]);
    adjpos2(i) = trapz([time(1:i)],[adjust2(1:i)]);
end
figure(3)
subplot(3,1,2)
hold on
plot(time,posG2)
plot(time,posE2)
plot(time,adjpos2)
title('Angular Position Vs Time, .7A .2Hz')
xlabel('Time (ms)')
ylabel('Radians')
legend('Gyro Output','Encoder Input','Adjusted Gyro Readings')
%third run testing mems gyro
Run3 = load('AM_Unit06_MEMS_COMP_RUN03');
Run3([1:31],:) = [];
time = Run3(:,1);
gyro_output = -Run3(:,2); %rads/s
input_rate = Run3(:,3); %rpm

input_rate = input_rate*(1/60)*(2*pi); %convert rpm to rads/s
figure(1)
subplot(3,1,3)
hold on
plot(input_rate,gyro_output)
title('Input rate VS Gyro output, .3A .7Hz')
xlabel('Encoder Input rad/ms')
ylabel('Gyro Output rad/ms')
legend('Data','Best fit line')

P = polyfit(input_rate,gyro_output,1);
y = polyval(P,input_rate);
bias = [bias P(2)];
sens = [sens P(1)];
plot(input_rate,y,'LineWidth',2)
legend('Data','Best fit line')

adjust3 = gyro_output*(1/sens(end)) - bias(end);
figure(2)
subplot(3,1,3)
hold on
plot(Run3(:,1),input_rate)
plot(Run3(:,1),gyro_output)
plot(Run3(:,1),adjust3)
title('Encoder Vs Gyro Vs Adjusted Angular Rates, .3A .7Hz')
xlabel('Time (ms)')
ylabel('Angular Velocity (rad/ms)')
legend('Input Rate','Gyro Output','Adjusted')
for i = 2:length(time)
    posG3(i) = trapz([time(1:i)],[gyro_output(1:i)]);
    posE3(i) = trapz([time(1:i)],[input_rate(1:i)]);
    adjpos3(i) = trapz([time(1:i)],[adjust3(1:i)]);
end
figure(3)
subplot(3,1,3)
hold on
plot(time,posG3)
plot(time,posE3)
plot(time,adjpos3)
title('Angular Position Vs Time, .3A .7Hz')
xlabel('Time (ms)')
ylabel('Radians')
legend('Gyro Output','Encoder Input','Adjusted Gyro Readings')

