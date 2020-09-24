%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 BALANCED AND UNBALANCED WHEEL LAB
% Date created:3/7/18
% Date modified:3/7/18
% 
% Author: Sage Herrin
% Takes time and angular velocity array and calculates and returns angular
% acceleration by take derivative of angular velocity (change in w/change
% in time).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha] = calc_alpha(time, theta_out, w_out)
%equal timesteps eliminates need for array of times, dt = 0.1 seconds =
%constant
alpha = [];
%define time step
dt = 0.1;
%takes difference between angular velocities and divides difference by time
%step
for i = 2:length(w_out)
    alpha(i - 1) = (w_out(i) - w_out(i - 1))/dt;
end
alpha = alpha';
% adjust vector length for plot compatability 
new_time = time(1:end - 1);
new_theta = theta_out(1:end - 1);
% plot data
plot(new_theta,alpha);
box off
grid minor
xlabel('Angular Position [rad]','FontSize',16,'Interpreter','latex')
ylabel('Angular Acceleration [rad/s$^2$]','FontSize',16,'Interpreter','latex')
title('unbalanced2.txt Data','FontSize',16,'Interpreter','latex');
