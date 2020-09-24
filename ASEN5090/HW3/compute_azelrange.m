function [az,el,range] = compute_azelrange(userECEF,satECEF)
LOS_ENU = compute_LOS_ENU(userECEF,satECEF);
az = atand(LOS_ENU(1)/LOS_ENU(2));
el = asind(LOS_ENU(3)/norm(LOS_ENU));
range = LOS_ENU(3);
% out = [az el range];
end