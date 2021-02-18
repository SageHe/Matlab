function Xs_ECI = stats_ECI(stats_ecef,t,theta0)
We = 7.2921158553e-5;
theta = theta0 + t*We; %Earth rotation angle in radians
C = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]; %DCM from ECEF to ECI

for i = 1:size(stats_ecef,1)
    Xs_ECI(i,1:3) = C*stats_ecef(i,:)';
%     Xs_ECI(i,4:6) = cross([0 0 We],Xs_ECI(i,1:3));
    Xs_ECI(i,4:6) = C*(cross([0 0 We],[stats_ecef(i,:)]))';


end
end