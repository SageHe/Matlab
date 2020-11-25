clear all; close all;clc
% read in given almanac data
[gps_alm,gps_alm_cell] = read_GPSyuma('YUMA310.ALM',2);
userpos_ecef = lla2ecef([40.01043021 -105.24388889 1614.976]); %negative 105 due to W not E

weeknum = 82; %week 82 without modulo 1024
tow = 591660;

for i = 1:size(gps_alm,1)
    [health(i),pos(i,:)] = broadcast_eph2pos(gps_alm,[weeknum tow],i);
    [az(i),el(i),range(i)] = compute_azelrange(userpos_ecef,pos(i,:));
end
