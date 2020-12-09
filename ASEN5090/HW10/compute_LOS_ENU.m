function out = compute_LOS_ENU(userECEF,satECEF)
vec = satECEF - userECEF;
lla_vec = ecef2lla(userECEF);
ECEF2ENU = calcECEF2ENU(lla_vec(1),lla_vec(2));
vec = ECEF2ENU*vec';
out = vec;
end