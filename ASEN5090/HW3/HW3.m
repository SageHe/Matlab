% ASEN 5090 Homework 3
% Author: Sage Herrin
% Date created: 9/22/20
clear all;close all;clc
%%Problem 1 -- Make a table of ECEF positions and latitude, longitude, and
%%height (wrt WGS84 ellipsoid) for all three user positions. 
NISTECEF = [-1288398 -4721697 4078625];
NISTLLA = ecef2lla(NISTECEF);

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
%% Problem 3 -- Write a function that takes as input ECEF opsitoin of a user on the ground and the ECEF
%position of a GPS satellite, and outputs the LOS unit vector pointing from
%the user to the satellite, expressed in ENU coords. wrt the user location
%on the Earth.
userpos_ECEF = [-1288398 -4721697 4078625];
satECEF = 2*userpos_ECEF;
LOS_ENU = compute_LOS_ENU(userpos_ECEF,satECEF);
LOS_ENU = LOS_ENU/norm(LOS_ENU);
%% Problem 4 -- Write a func. that coputes the geometrical range (m), azimuth (deg), and elevation (deg) of a 
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
userpos_ECEF = [-1288398 -4721697 4078625]; %Put into km from m
for i = 1:size(Data,1)
    [AZ(i), EL(i), RANGE(i)] = compute_azelrange(userpos_ECEF, [Data(i,2:4)]*1000);    
end
vis = EL>0;
PRN = Data(vis,1);
AZ = AZ(vis);
EL = EL(vis);
plotAzEl(AZ,EL,PRN)