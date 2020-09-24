%Complab 4
%Author: Sage Herrin
%Created: 11/9/18
%Function that solves PLLT for thick airfoils, evaluates left wing,
%-b/2 < y <0
%Enter geometric AOA in degrees, convert to radians later 
function [e, c_L, c_Di] = PLLT[b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%Acquire necessay alpha values from lab 3 and given values 
%Use given geometric AOA values-degrees & map to a(theta)
geo_t = geo_t*(pi/180);
geo_r = geo_r*(pi/180);
%Get zero_lift AOA from lab 3
% aero_t = .0017; %rads - technically should be 0, try tweaking later
% aero_r = -2.1184; %rads
%Get cross-sectional lift slopes from lab 3
% a0_t = 7.0459;
% a0_r = 7.0326;
%Calculate various AOAs and map to theta
for i = 1:N
    theta(i) = (i*pi)/(2*N);
end
y = (-b/2)*cos(theta);
geo_AOA = geo_r - (geo_t - geo_r)*y/(b/2);
AOA_0L = aero_r - (aero_t - aero_r)*y/(b/2);
end