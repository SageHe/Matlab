function Xspos_ECEF = stats_ECEF(stats_ECI,t,theta0)
We = 7.2921158553e-5;
theta = theta0 + t*We; %Earth rotation angle in radians
C = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]; %DCM from ECI to ECEF

for i = 1:size(stats_ECI,1)
    Xspos_ECEF(i,1:3) = C'*stats_ECI(i,:)';
end
end