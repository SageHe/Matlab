% ASEN 5090 Homework 3
% Author: Sage Herrin
% Date created: 9/22/20
clear all;close all;clc
%%Problem 1 -- Make a table of ECEF positions and latitude, longitude, and
%%height (wrt WGS84 ellipsoid) for all three user positions. 
NISTECEF = [-1288398 -4721697 4078625];
NISTLLA = ecef2lla(NISTECEF);
EQUAECEF = lla2ecef([0 NISTLLA(2) 0]);
ENNOECEF = lla2ecef([89 NISTLLA(2) 0]);
%% Problem 2 -- Write a function that takes as input a reference latitude and longitude in degrees 
%and returns the transformation matrix to rotate from ECEF to ENU

%Rotation matrix from ECEF to ENU
lambda = NISTLLA(2);
phi = NISTLLA(1);
ECEF2ENU = calcECEF2ENU(lambda,phi);
% ECEF2ENU = [-sind(lambda) cosd(lambda) 0;...
%             -sind(phi)*cosd(lambda) -sind(phi)*sind(lambda) cosd(phi);...
%             cosd(phi)*cosd(lambda) cosd(phi)*sind(lambda) sind(phi)];
NISTENU = ECEF2ENU*NISTECEF'; 
%% Problem 3 -- Write a function that takes as input ECEF positoin of a user on the ground and the ECEF
%position of a GPS satellite, and outputs the LOS unit vector pointing from
%the user to the satellite, expressed in ENU coords. wrt the user location
%on the Earth.
userpos_ECEF = [-1288398 -4721697 4078625];
satECEF = 2*userpos_ECEF;
LOS_ENU = compute_LOS_ENU(userpos_ECEF,satECEF);
LOS_ENU = LOS_ENU/norm(LOS_ENU);
%% Problem 4 -- Write a func. that computes the geometrical range (m), azimuth (deg), and elevation (deg) of a 
% GPS satellite based on the ECEF locations of a user and GPS satellite.
% Make up an example to check your func.
[AZ, EL, RANGE] = compute_azelrange(userpos_ECEF, satECEF);
%% Problem 5 -- Find the GPS satellite visibility based on satellite and user positions
Data = load('gpspos20200901.txt');
% loop through each user location and:
    %{
    i. Compute az, el, and range for each GPS sat. for user pos.
    ii. Determine which sats. are visible
    iii. Make a table with PRN, AZ, EL, Range for visible. sats.
    iv. Plot visible sats. on skyplot: plotazel(AZ,EL,PRN)
    %}
% For first user location (NIST)
clear AZ EL RANGE
userpos = [NISTECEF;EQUAECEF;ENNOECEF]; %[-1288398 -4721697 4078625]; %Put into km from m
for j = 1:3
    userpos_ECEF = userpos(j,:);
    for i = 1:size(Data,1)
        [AZ(i), EL(i), RANGE(i)] = compute_azelrange(userpos_ECEF, [Data(i,2:4)]*1000);    
    end
    vis = EL>0;
    PRN = Data(vis,1);
    AZ = AZ(vis);
    EL = EL(vis);
    figure
    plotAzEl(AZ,EL,PRN)
end
%% Problem 6 -- Write a script to compute the GPS sat. positions for the entire day using the GPS almanac data. Script should:
    %{
        a. Use function read_PGSyuma to load the provided almanac for sept.
        1 2020
        b. Set up a time of week (in secs) vector that starts at the
        beginning of 9/1/2020 and runs for 24 hours at reasonable spacing
        (1 min)
        c. Compute the positions of the GPS sats. based on the almanac
        using the function broadcast2pos provided
    %}
[gps_ephem,gps_ephem_cell] = read_GPSyuma('YUMA245.ALM',2);
% set up time of week vector in seconds
tvec = gps_ephem_cell{1,1}.Toe:60:gps_ephem_cell{1,1}.Toe+24*3600;
% compute positions of GPS satellites based on almanac using function
% broadcast2pos
for i = 1:32
    [health(:,i),pos(:,:,i)] = broadcast2pos(gps_ephem,[73*ones(length(tvec),1) tvec'],i);
end
%% Problem 7 -- Write a script to compute and plot the satellite visibility over the entire data, for each of the 3
% ground locations and the satellite positions predicted using the almanac
% data
% convert gps satellite pos. data into az el range
for h = 1:3
    for i = 1:size(pos,3)
        for j = 1:size(pos,1)
            [az(j,i,h) el(j,i,h) range(j,i,h)] = compute_azelrange(userpos(h,:),pos(j,:,i));
        end
    end
end
% for i = 1:size(pos,1)
% note: not very efficient, not passing vectors into plotAzEl func       
for h = 1:3
figure
    for i = 1:size(pos,3)
        for j = 1:size(pos,1)
            if el(j,i,h)>0
            plotAzEl(az(j,i,h),el(j,i,h),zeros(1,size(pos,3)))
            end
            [j i]
        end
    end
end