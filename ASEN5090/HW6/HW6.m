% Homework 6: Least Squares Estiamte for GPS Position
% Author: Sage Herrin
% Last Modified: 10/24/20
%Housekeeping
clear all;close all;clc
%Part 1 -- Use given NIST location as the initial guess of the user position
NISTECEF = [-1288398.360 -4721697.040 4078625.5];
C =  2.99792458e8; %m/s
%part 2 -- Read first epoch of observation records. For each of the observed sats., print and check the following values:
    %{
        a. iono free pseudorange observable PIF
        b. expected range in meters
        c. elevation and azimuth in degrees
        d. sat. clock correction in meters
        e. relativistic correction in meters
        f. simple tropo. model in meters
    %}
ephem_data = read_clean_GPSbroadcast('brdc2450.20n',true);

% rinex_data = read_rinex_obs8('nist2450.20o');
load('epoch1');
epoch1 = data;
rinex_tvec = epoch1(1,2);
Weeknum_vec = epoch1(1);
% epoch1 = rinex_data.data(1:10,:);
f1 = 1575.42e6;
f2 = 1227.6e6;
for i = 1:size(epoch1,1)
    [PIF(i),~] = ionocorr(epoch1(i,4),f1,epoch1(i,8),f2);
    [health_eph(i,:),pos_eph(i,:)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_tvec],epoch1(i,3));
    [az,el,range] = compute_azelrange(NISTECEF,pos_eph(i,:));
    for j = 1:2
    Tt = rinex_tvec' - (range./C);
    % Tt = pos2_eph - (range2'./C);
    %Compute satellite position at Tt in ECEF at Tt based on broadcast
    %ephemeris
    [health_ephTt,pos_ephTt(i,:),bsv(i),relsv(i)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec Tt'],epoch1(i,3));
    %Rotation rate of Earth
    We = 7.2921151467e-5; %rads/sec, (hypertextbook)
    phi = We*(rinex_tvec' - Tt);
    for k = 1:length(Tt)
        ECEF_rot(i,:) = [cos(phi(k)) sin(phi(k)) 0;...
                        -sin(phi(k)) cos(phi(k)) 0;...
                        0 0 1]*pos_ephTt(i,:)';
    R(k) = norm(ECEF_rot(i,:) - NISTECEF);               
    end
    diff = abs(R - range);
    range = R;
    end
    exp_range(i) = range;
    [az_corr(i),el_corr(i),~] = compute_azelrange(NISTECEF,ECEF_rot(i,:));
    tropo(i) = tropomodel(el_corr(i),2);
end
%% Part 3 -- construct and print A or G matrix using ECEF coords of GPS sats
% and ovserver
for i = 1:10
    G(i,1) = -((ECEF_rot(i,1)-NISTECEF(1))/exp_range(i));
    G(i,2) = -((ECEF_rot(i,2)-NISTECEF(2))/exp_range(i));
    G(i,3) = -((ECEF_rot(i,3)-NISTECEF(3))/exp_range(i));
    G(i,4) = 1;
end
%% Part 4 -- Compute and print prefit residuals correcting for ionospheric dealy, sat. clock, relativity, and tropo. delay
dy = PIF - (exp_range - bsv - relsv + tropo);
dy = dy';
%% Part 5 -- plot residuals as a function of the sat. elevation angle
figure
plot(el_corr,dy,'*')
ylim([0 100])
%% Part 6 -- Form least squares solution for delta x state vector
dx = (G'*G)\G'*dy;
%% Part 7 -- Use initial estimate that is way off of correct user position.
% Compute A matix and pre-fit residuals based on this assumed pos.. Compute
% least squares correction and iterate soln. 5 times. 
userECEF = [0 0 0];
b = 0;
db = 0;
for h = 1:5
    for i = 1:size(epoch1,1)
        [PIF(i),~] = ionocorr(epoch1(i,4),f1,epoch1(i,8),f2);
        [health_eph(i,:),pos_eph(i,:)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_tvec],epoch1(i,3));
        [az,el,range] = compute_azelrange(userECEF,pos_eph(i,:));
        for j = 1:2
        Tt = rinex_tvec' - (range./C);
        % Tt = pos2_eph - (range2'./C);
        %Compute satellite position at Tt in ECEF at Tt based on broadcast
        %ephemeris
        [health_ephTt,pos_ephTt(i,:),bsv(i),relsv(i)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec Tt'],epoch1(i,3));
        %Rotation rate of Earth
        We = 7.2921151467e-5; %rads/sec, (hypertextbook)
        phi = We*(rinex_tvec' - Tt);
        for k = 1:length(Tt)
            ECEF_rot(i,:) = [cos(phi(k)) sin(phi(k)) 0;...
                            -sin(phi(k)) cos(phi(k)) 0;...
                            0 0 1]*pos_ephTt(i,:)';
        R(k) = norm(ECEF_rot(i,:) - userECEF);               
        end
        diff = abs(R - range);
        range = R;
        end
        exp_range(i) = range;
        [az_corr(i),el_corr(i),~] = compute_azelrange(userECEF,ECEF_rot(i,:));
        tropo(i) = tropomodel(el_corr(i),2);
    end

    for i = 1:10
        G(i,1) = -((ECEF_rot(i,1)-userECEF(1))/exp_range(i));
        G(i,2) = -((ECEF_rot(i,2)-userECEF(2))/exp_range(i));
        G(i,3) = -((ECEF_rot(i,3)-userECEF(3))/exp_range(i));
        G(i,4) = 1;
    end
    b = b + db;
    dy = PIF - (exp_range + b - bsv - relsv + tropo);
    dy = dy';

    dx = (G'*G)\G'*dy;
    userECEF = userECEF + dx(1:3);
    db = dx(4);
end