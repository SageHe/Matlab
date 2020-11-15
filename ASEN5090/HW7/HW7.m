clear all;close all;clc
%% Models and preparatoin for solution
% Set a priori locations for NIST and USN8
NISTECEF = [-1288398.360 -4721697.040 4078625.500];
USN8ECEF = [1112161.9919 -4842855.0890  3985497.3152];
C =  2.99792458e8; % speed of light,m/s
% Compute transformation matrix from ecef to enu using provided locations
% NISTLLA = ecef2lla(NISTECEF);
% USN8LLA = ecef2lla(USN8ECEF);
% C_ecef2enu_nist = calcECEF2ENU(NISTLLA(1),NISTLLA(2));
% C_ecef2enu_usn8 = calcECEF2ENUUSNA8LLA(1),USN8LLA(2));
% Set up code to loop through entire obs file, start with first few epochs while debugging
rinex_data = read_rinex_obs8('nist2450.20o',[1:32]);
ephem_data = read_clean_GPSbroadcast('brdc2450.20n',true);
Weeknum_vec = rinex_data.data(1,1);
epoch = 1;
satcnt = 1;
data = [];
for i = 1:size(rinex_data.data,1)-1
    if rinex_data.data(i,2) == rinex_data.data(i+1,2)
        %Calculate all the stuff
        [health_eph,pos_eph,bsv,relsv] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_data.data(i,2)],rinex_data.data(i,3));
        rinex_tvec = rinex_data.data(i,2);
%         data(satcnt,1:3,epoch) = pos_eph;
%         data(satcnt,4,epoch) = bsv;
%         data(satcnt,5,epoch) = relsv;
        f1 = 1575.42e6;
        f2 = 1227.6e6;
%         for i = 1:size(epoch1,1)
        % calc expected range
        [PIF,~] = ionocorr(rinex_data.data(i,4),f1,rinex_data.data(i,8),f2);
%         [health_eph(i,:),pos_eph(i,:)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_tvec],epoch1(i,3));
        [az,el,range] = compute_azelrange(NISTECEF,pos_eph);
        if el < 10
            data(satcnt,:,epoch) = 0;
            satcnt = satcnt+1;
            continue
        end
        for j = 1:2
        Tt = rinex_tvec' - (range./C);
        % Tt = pos2_eph - (range2'./C);
        %Compute satellite position at Tt in ECEF at Tt based on broadcast
        %ephemeris
        [health_ephTt,pos_ephTt,bsv,relsv] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec Tt'],rinex_data.data(i,3));
        data(satcnt,1:3,epoch) = pos_ephTt;
        data(satcnt,4,epoch) = bsv;
        data(satcnt,5,epoch) = relsv;
        %Rotation rate of Earth
        We = 7.2921151467e-5; %rads/sec, (hypertextbook)
        phi = We*(rinex_tvec' - Tt);
        for k = 1:length(Tt)
            ECEF_rot(i,:) = [cos(phi(k)) sin(phi(k)) 0;...
                            -sin(phi(k)) cos(phi(k)) 0;...
                            0 0 1]*pos_ephTt';
        R(k) = norm(ECEF_rot(i,:) - NISTECEF);               
        end
        diff = abs(R - range);
        range = R;
        end
%         exp_range(i) = range;
        [az_corr(i),el_corr(i),~] = compute_azelrange(NISTECEF,ECEF_rot(i,:));
        tropo = tropomodel(el_corr(i),2);
        data(satcnt,6,epoch) = range; %expected range of sat.
        data(satcnt,7,epoch) = tropo;
        data(satcnt,8,epoch) = PIF - (range - bsv - relsv + tropo);
        data(satcnt,9,epoch) = el_corr(i);
        satcnt = satcnt+1;
%          end
    else
        %Calculate all the stuff then increment epoch counter
        [health_eph,pos_eph] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_data.data(i,2)],rinex_data.data(i,3));
        rinex_tvec = rinex_data.data(i,2);
        data(satcnt,1:3,epoch) = pos_eph;
        data(satcnt,4,epoch) = bsv;
        data(satcnt,5,epoch) = relsv;
        f1 = 1575.42e6;
        f2 = 1227.6e6;
%         for i = 1:size(epoch1,1)
        % calc expected range
        [PIF,~] = ionocorr(rinex_data.data(i,4),f1,rinex_data.data(i,8),f2);
%         [health_eph(i,:),pos_eph(i,:)] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec rinex_tvec],epoch1(i,3));
        [az,el,range] = compute_azelrange(NISTECEF,pos_eph);
        if el < 10
            data(satcnt,:,epoch) = 0;
            satcnt = 1;
            epoch = epoch+1;
            continue
        end
        for j = 1:2
        Tt = rinex_tvec' - (range./C);
        % Tt = pos2_eph - (range2'./C);
        %Compute satellite position at Tt in ECEF at Tt based on broadcast
        %ephemeris
        [health_ephTt,pos_ephTt,bsv,relsv] = broadcast_eph2pos_etc(ephem_data,[Weeknum_vec Tt'],rinex_data.data(i,3));
        data(satcnt,1:3,epoch) = pos_eph;
        data(satcnt,4,epoch) = bsv;
        data(satcnt,5,epoch) = relsv;
        %Rotation rate of Earth
        We = 7.2921151467e-5; %rads/sec, (hypertextbook)
        phi = We*(rinex_tvec' - Tt);
        for k = 1:length(Tt)
            ECEF_rot(i,:) = [cos(phi(k)) sin(phi(k)) 0;...
                            -sin(phi(k)) cos(phi(k)) 0;...
                            0 0 1]*pos_ephTt';
        R(k) = norm(ECEF_rot(i,:) - NISTECEF);               
        end
        diff = abs(R - range);
        range = R;
        end
%         exp_range(i) = range;
        [az_corr(i),el_corr(i),~] = compute_azelrange(NISTECEF,ECEF_rot(i,:));
        tropo = tropomodel(el_corr(i),2);
        data(satcnt,6,epoch) = range; %expected range of sat.
        data(satcnt,7,epoch) = tropo; %tropo delay
        data(satcnt,8,epoch) = PIF - (range - bsv - relsv + tropo);
        data(satcnt,9,epoch) = el_corr(i);
        satcnt = 1;
        epoch = epoch+1;
%         end
    end
end
%% Least Squares Point Solutions
for i = 1:size(data,3)
    temp = data(:,:,i);
    for j = 1:size(data,1)
        if isnan(temp(j,:))
            temp(j,:) = [];
            continue
        end
        if ~any(temp(j,:))
            temp(j,:) = [];
            continue
        end
        G(i,1) = -((ECEF_rot(i,1)-NISTECEF(1))/exp_range(i));
        G(i,2) = -((ECEF_rot(i,2)-NISTECEF(2))/exp_range(i));
        G(i,3) = -((ECEF_rot(i,3)-NISTECEF(3))/exp_range(i));
        G(i,4) = 1;
            
            
            
            
            