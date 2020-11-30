close all;clc
%% Part 1 -- Visibility Prediction
% read in given almanac data
data = load('HW9data.mat');
[gps_alm,gps_alm_cell] = read_GPSyuma('YUMA310.ALM',2);
userpos_ecef = lla2ecef([40.01043021 -105.24388889 1614.976]); %negative 105 due to W not E

weeknum = 2130; %week 82 without modulo 1024
tow = 591660;
prn = [1:14 16:32];
for i = prn
    [health(i),pos(i,:)] = broadcast_eph2pos(gps_alm,[weeknum tow],i);
    [az(i),el(i),range(i)] = compute_azelrange(userpos_ecef,pos(i,:));
end
%% Part 2 -- Carrier Wipeoff/Code Correlation
Sr = data.gpsdata;

tvec = linspace(0,1e-3,5001);
prn31 = CA_Gen(31,1024);
% prn31 = flip(prn31);
prn31 = (prn31==0)*1 + (prn31==1)*(-1);
chiprate = 1.023e6;
Tchip = 1/chiprate;
ind = floor(tvec/Tchip)+1;
sampprn = prn31(ind);

IF  = 1.25e6; %intermediate freq.
fd = 0; %Doppler fequency 

theta = 2*pi*(IF+fd)*tvec;
delay = 9;
S = comp_corr(Sr,delay,theta,numel(tvec),sampprn);
%% Part 3 -- Delay/Doppler Grid Search
delayax = 0:1:5000; %change this to 5001 and see if it still works
doppax = -10e3:500:10e3;

%% Part 4
% found 32, 1, 22, 10