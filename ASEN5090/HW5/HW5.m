%Author: Sage Herrin
%Created: 10/15/20
clear all;close all;clc
%Load in ephem data
ephem_data = read_clean_GPSbroadcast('brdc2450.20n',true);
%Load in rinex file
rinex_data = read_rinex_obs8('nist2450.20o',1);
P2_vals = rinex_data.data(:,8);
rinex_tvec = rinex_data.data(:,2);

Weeknum = ephem_data(1,19);
dow = 2;
Tow = (0:30:24*3600) + dow*24*3600;
Tow_hr = (Tow - dow*24*3600)/3600;

Weeknum_vec = Weeknum*ones(size(rinex_tvec,1),1);
[health_eph,pos_eph] = broadcast_eph2pos(ephem_data,[Weeknum_vec rinex_tvec],1);
NISTECEF = [-1288398.360 -4721697.040 4078625.5];

for i = 1:length(rinex_tvec)
    [az(i),el(i),range(i)] = compute_azelrange(NISTECEF,pos_eph(i,:));
end

C =  2.99792458e8; %m/s

orig_range = range;

for j = 1:3
    Tt = rinex_tvec' - (range./C);
    % Tt = pos2_eph - (range2'./C);
    %Compute satellite position at Tt in ECEF at Tt based on broadcast
    %ephemeris
    [health_ephTt,pos_ephTt,bsv,relsv] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec Tt'],1);
    %Rotation rate of Earth
    We = 7.2921151467e-5; %rads/sec, (hypertextbook)
    phi = We*(rinex_tvec' - Tt);
    for i = 1:length(Tt)
        ECEF_rot(i,:) = [cos(phi(i)) sin(phi(i)) 0;...
                        -sin(phi(i)) cos(phi(i)) 0;...
                        0 0 1]*pos_ephTt(i,:)';
    R(i) = norm(ECEF_rot(i,:) - NISTECEF);               
    end
    diff = abs(R - range);
    range = R;
end
dPR0 = P2_vals - range';
rinex_tvec = (rinex_tvec - dow*86400)/3600;
figure
plot(rinex_tvec,dPR0)
grid on
figure
plot(rinex_tvec,bsv)
grid on
dPR1 = P2_vals - (range' - bsv');
figure
plot(rinex_tvec,dPR1)
grid on
dPR2 = P2_vals - (range' - bsv' - relsv');
figure
plot(rinex_tvec,relsv)
grid on
figure
plot(rinex_tvec,dPR2)
grid on 
tropo = tropomodel(el,2);
figure
plot(rinex_tvec,tropo)
grid on
dPR3 = P2_vals - (range' - bsv' - relsv' + tropo');
figure
plot(rinex_tvec,dPR3)
grid on
f1 = 1575.42;
f2 = 1227.6;
for i = 1:size(rinex_data.data(:,4))
    C1(i) = rinex_data.data(i,4);
    P2(i) = rinex_data.data(i,8);
    [PRIF(i),iono(i)] = ionocorr(C1(i),f1,P2(i),f2);
end
figure
plot(rinex_tvec,iono)
grid on
dPR4 = PRIF' - (range' - bsv' - relsv' + tropo');
figure
plot(rinex_tvec,dPR4)
grid on