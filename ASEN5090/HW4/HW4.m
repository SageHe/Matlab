%% ASEN 5090 Homework 4 -- Computing Expected Range
% Author :Sage Herrin, with help and input from Lara Lufkin
% Date Created: 10/1/20
close all;clc
%% Part 1 -- Determine name of full boradcase ephemeris file for Sept. 1, 2020
filename = "brdc2450.20n";
%% Part 2 -- load ephemeris data into numerical array
ephem_data = read_clean_GPSbroadcast('brdc2450.20n',true);
%% Part 3 -- Compute pos. of a GPS sat. based on the ephemeris data, for specified set of times.
Weeknum = ephem_data(1,19);
dow = 2;
Tow = (0:30:24*3600) + dow*24*3600;
Tow_hr = (Tow - dow*24*3600)/3600;
Weeknum_vec = Weeknum*ones(size(Tow,2),1);
[health,pos] = broadcast_eph2pos(ephem_data,[Weeknum_vec Tow'],19);
%% Part 4 -- plot ECEF coords computed based on almanac for two sats. and compare to to ECEF coords. computed using new eph2pos function 
%for same times.
%first input almanac data and create time vector
[gps_alm,gps_alm_cell] = read_GPSyuma('YUMA245.ALM',2);
% tvec_alm = 0:30:24*3600;
% tvec_alm = gps_ephem_cell{1,1}.Toe:30:gps_ephem_cell{1,1}.Toe+24*3600;
Weeknumvec_alm = Weeknum*ones(size(Tow,2),1);
[health2_alm,pos2_alm] = broadcast2pos(gps_alm,[Weeknumvec_alm Tow'],2);
[health24_alm,pos24_alm] = broadcast2pos(gps_alm,[Weeknumvec_alm Tow'],24);

[health2_eph,pos2_eph] = broadcast_eph2pos(ephem_data,[Weeknum_vec Tow'],2);
[health24_eph,pos24_eph] = broadcast_eph2pos(ephem_data,[Weeknum_vec Tow'],24);

figure 
hold on
plot(Tow_hr',pos2_alm)
plot(Tow_hr',pos2_eph,'--')
xlabel('Time (hrs)')
ylabel('ECEF Position (m)')
title('ECEF Position VS Time, PRN 2')
figure
hold on
plot(Tow_hr',pos24_alm)
plot(Tow_hr',pos24_eph,'--')
xlabel('Time (hrs)')
ylabel('ECEF Position (m)')
title('ECEF Position VS Time, PRN 24')

satdiff2 = pos2_alm - pos2_eph;
satdiff24 = pos24_alm - pos24_eph;
figure
plot(Tow_hr',satdiff2)
xlabel('Time (hrs)')
ylabel('Position Difference (m)')
title('Difference in Almanac VS Ephemeris Position, PRN 2')
figure
plot(Tow_hr',satdiff24)
xlabel('Time (hrs)')
ylabel('Position Difference (m)')
title('Difference in Almanac VS Ephemeris Position, PRN 24')
%% Part 5 -- find the range, az., and el. from sats. to the NIST IGS site for entire day. Use GPS positions calculated using new ephemeris function.
NISTECEF = [-1288398.360 -4721697.040 4078625.5];
for i = 1:length(Tow)
    [az2(i),el2(i),range2(i)] = compute_azelrange(NISTECEF,pos2_eph(i,:));
    [az24(i),el24(i),range24(i)] = compute_azelrange(NISTECEF,pos24_eph(i,:));
end
figure
hold on
plot(Tow_hr,az2)
plot(Tow_hr,az24)
xlabel('Time (hrs)')
ylabel('Azimuth (degrees)')
title('Azimuth Angle VS Time')
legend('PRN 2','PRN 24')
figure 
hold on 
plot(Tow_hr,el2)
plot(Tow_hr,el24)
xlabel('Time (hrs)')
ylabel('Elevation (degrees)')
title('Elevation Angle VS Time')
legend('PRN 2','PRN 24')
figure 
hold on
plot(Tow_hr,range2)
plot(Tow_hr,range24)
xlabel('Time (hrs)')
ylabel('Range (m)')
title('Range VS Time')
legend('PRN 2','PRN 24')
%% Part 6 -- Recompute the range accounting for signal travel time and coord. frame rotation as described in attachment.
C =  2.99792458e8; %m/s
%Compute time of transmission
orig_range2 = range2;
for j = 1:3
    Tt = Tow - (range2./C);
    % Tt = pos2_eph - (range2'./C);
    %Compute satellite position at Tt in ECEF at Tt based on broadcast
    %ephemeris
    [health2_ephTt,pos2_ephTt] = broadcast_eph2pos(ephem_data,[Weeknum_vec Tt'],2);
    %Rotation rate of Earth
    We = 7.2921151467e-5; %rads/sec, (hypertextbook)
    phi = We*(Tow - Tt);
    for i = 1:length(Tt)
        ECEF_rot(i,:) = [cos(phi(i)) sin(phi(i)) 0;...
                        -sin(phi(i)) cos(phi(i)) 0;...
                        0 0 1]*pos2_ephTt(i,:)';
    R(i) = norm(ECEF_rot(i,:) - NISTECEF);               
    end
    diff = abs(R - range2);
    range2 = R;
end
plot(Tow_hr,R)
rangediff = range2 - orig_range2;
figure
hold on
plot(Tow_hr,rangediff)
%% Part 7 -- Extract all measured C1 and P2 pseudoranges for satellite that is chosen. Plot pseudoranges vs time in hours of the day
rinex_data = read_rinex_obs8('nist2450.20o',1);
C1_2_vals = rinex_data.data(:,4);
P2_2_vals = rinex_data.data(:,8);
rinex_tvec_2 = rinex_data.data(:,2);
% C1_24_vals = rinex_data.data(find(rinex_data.data(:,3) == 24),4);
% P2_24_vals = rinex_data.data(find(rinex_data.data(:,3) == 24),8);
% rinex_tvec_24 = rinex_data.data(find(rinex_data.data(:,3) == 24),2);

figure
hold on
plot(rinex_tvec_2,C1_2_vals)
plot(rinex_tvec_2,P2_2_vals)
% figure
% hold on
% plot(rinex_tvec_24,C1_24_vals)
% plot(rinex_tvec_24,P2_24_vals)
%% Part 8 -- Compare measured pseudoranges from rinex file to predicted preudoranges
clear range2 R diff health2_ephTt pos2_ephTt pos2_eph ECEF_rot
Weeknum_vec = Weeknum*ones(size(rinex_tvec_2,1),1);
[health2_eph,pos2_eph] = broadcast_eph2pos(ephem_data,[Weeknum_vec rinex_tvec_2],1);

for i = 1:length(rinex_tvec_2)
    [az2(i),el2(i),range2(i)] = compute_azelrange(NISTECEF,pos2_eph(i,:));
end

orig_range2 = range2;

for j = 1:3
    Tt = rinex_tvec_2' - (range2./C);
    % Tt = pos2_eph - (range2'./C);
    %Compute satellite position at Tt in ECEF at Tt based on broadcast
    %ephemeris
    [health2_ephTt,pos2_ephTt] = broadcast_eph2pos(ephem_data,[Weeknum_vec Tt'],1);
    %Rotation rate of Earth
    We = 7.2921151467e-5; %rads/sec, (hypertextbook)
    phi = We*(rinex_tvec_2' - Tt);
    for i = 1:length(Tt)
        ECEF_rot(i,:) = [cos(phi(i)) sin(phi(i)) 0;...
                        -sin(phi(i)) cos(phi(i)) 0;...
                        0 0 1]*pos2_ephTt(i,:)';
    R(i) = norm(ECEF_rot(i,:) - NISTECEF);               
    end
    diff = abs(R - range2);
    range2 = R;
end
figure
hold on
diffrange_2 = P2_2_vals - range2';
plot(rinex_tvec_2,P2_2_vals)
plot(rinex_tvec_2,range2)
figure
plot(rinex_tvec_2,diffrange_2)