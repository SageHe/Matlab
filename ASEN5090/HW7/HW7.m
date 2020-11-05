% Set a priori locations for NIST and USN8
NISTECEF = [-1288398.360 -4721697.040 4078625.500];
USN8ECEF = [1112161.9919 -4842855.0890  3985497.3152];
% Compute transformation matrix from ecef to enu using provided locations
NISTLLA = ecef2lla(NISTECEF);
%USN8LLA = ecef2lla(USN8ECEF);
C_ecef2enu = calcECEF2ENU(NISTLLA(1),NISTLLA(2));
% Set up code to loop through entire obs file, start with first few epochs while debugging
rinex_data = read_rinex_obs8('USN82450.20o',[1:32],10);
ephem_data = read_clean_GPSbroadcast('brdc2450.20n',true);
