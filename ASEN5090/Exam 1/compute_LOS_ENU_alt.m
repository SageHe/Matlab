function out = compute_LOS_ENU_alt(userECEF,satECEF)
vec = satECEF - userECEF;
lla_vec = ecef2lla(userECEF,0,6380000);
ECEF2ENU = calcECEF2ENU(lla_vec(1),lla_vec(2));
vec = ECEF2ENU*vec';
out = vec;
end