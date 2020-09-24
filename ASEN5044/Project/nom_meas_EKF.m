function [Y] = nom_meas_EKF(X,X_gt,t,j)
Re = 6378;
we = (2*pi)/86400;
theta0 = (j-1)*(pi/6);
Xi = Re*cos(we*t+theta0);
Yi = Re*sin(we*t+theta0);
Xidot = -Re*sin(we*t+theta0)*we;
Yidot = Re*cos(we*t+theta0)*we;
theta_i = atan2(Yi,Xi);
rho = sqrt((X(1)-Xi)^2+(X(3)-Yi)^2);
rhodot = ((X(1)-Xi)*(X(2)-Xidot)+(X(3)-Yi)*(X(4)-Yidot))/rho;
phi = atan2((X(3)-Yi),(X(1)-Xi));
phi_check = atan2((X_gt(3)-Yi),X_gt(1)-Xi);
lower1 = -pi/2+theta_i;
upper1 = pi/2 +theta_i;
lower2 = lower1;
upper2 = upper1;
if -pi/2+theta_i < -pi
    lower1 = -pi/2+theta_i+2*pi;
    upper1 = pi;
    lower2 = -pi;
    upper2 = pi/2+theta_i;
else
if  pi/2+theta_i > pi
    lower1 = -pi;
    upper1 = pi/2+theta_i-2*pi;
    lower2 = -pi/2+theta_i;
    upper2 = pi;
end
end
if (phi_check >= lower1 && phi_check <= upper1)|| (phi_check >= lower2 && phi_check <= upper2)
    Y = [rho;rhodot;phi];
else
    Y = NaN(3,1);
end
% if phi < ((-pi/2)+theta_i) || phi > ((pi/2)+theta_i)
%     Y = nan(3,1);
% else
%     Y = [rho;rhodot;phi];
% end
end