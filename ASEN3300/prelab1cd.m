%% good housekeeping :)
close all
clear all
clc

%% establish constants
minV=0;
maxV=3.3;
voltages=[0 .25 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25];

%% 4 Bits
bits=4;

[outbin, binout]=Voltage2Bin(minV, maxV, bits, voltages);

figure
plot(voltages,outbin,'s')
title('Bins per Voltage 4 Bits')
ylabel('Bin Number')
xlabel('Voltage Value [V]')

%% 8 Bits

bits=8;

[outbin, binout]=Voltage2Bin(minV, maxV, bits, voltages);

figure
plot(voltages,outbin,'s')
title('Bins per Voltage 8 Bits')
ylabel('Bin Number')
xlabel('Voltage Value [V]')

%% 12 Bits

bits=12;

[outbin, binout]=Voltage2Bin(minV, maxV, bits, voltages);

figure
plot(voltages,outbin,'s')
title('Bins per Voltage 12 Bits')
ylabel('Bin Number')
xlabel('Voltage Value [V]')


%% Part d
t=linspace(0,2*pi,100);
v=1.65+3.3/2.*sin(t);
bits=12;
minV=1.65-3.3/2;
maxV=1.65+3.3/2;
[outbin, binout]=Voltage2Bin(minV, maxV, bits, v);


figure
plot(1:length(outbin),outbin)

title('Bin Number vs Array number Sinusoidal Input')
ylabel('Bin Number')
xlabel('Array Number')