function y = predictmeas(x_sc,x_obs,t)
R = [x_sc(1) x_sc(2) x_sc(3)];
% Rs = [x_obs(1) x_obs(2) x_obs(3)];
V = [x_sc(4) x_sc(5) x_sc(6)];
% Vs = [x_obs(4) x_obs(5) x_obs(6)];
y = zeros(1,2*size(x_obs,1));
%Transform station position and velocity into ECI frame for valid
%calculations
% We = 7.2921158553e-5;
% theta0 = deg2rad(122);
% theta = theta0 + We*t;
% C = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]; %DCM from ECEF to ECI
% 
% Rs_ECI = C*Rs;
% Vs_ECI = cross([0 0 We],Rs_ECI);
for i = 1:size(x_obs,1)
    rho = norm(R - x_obs(i,1:3));
    y(2*i-1) = rho;
    y(2*i) = (dot((R - x_obs(i,1:3)),(V - x_obs(i,4:6))))/rho;
    % rho = norm(R - Rs);
    % rhodot = (dot((R - Rs),(V - Vs))/rho);
end
