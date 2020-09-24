%Complab 4
%Author: Sage Herrin
%Created: 11/9/18
%Function that solves PLLT for thick airfoils, evaluates left wing,
%-b/2 < y <0
%Enter geometric AOA in degrees, convert to radians later 
function [e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%Acquire necessay alpha values from lab 3 and given values 
%Use given geometric AOA values-degrees & map to a(theta)
% geo_t = geo_t*(pi/180);
% geo_r = geo_r*(pi/180);
% aero_t = aero_t*(pi/180);
% aero_r = aero_r*(pi/180);

%Get zero_lift AOA from lab 3
% aero_t = .0017 deg, ~0 rads - technically should be 0, try tweaking later
% aero_r = -2.1184 deg, -.0370 rads
%Get cross-sectional lift slopes from lab 3
% a0_t = 7.0459;
% a0_r = 7.0326;
%Calculate various AOAs and map to theta
% theta(1) = 0;
for i = 1:N
    theta(i) = (i*pi)/(2*N);
end
% theta = [pi/3 pi/2];
y = (-b/2)*cos(theta);
geo_AOA = geo_r - (geo_r - geo_t).*cos(theta); %y/(b/2);
AOA_0L = aero_r - (aero_r - aero_t).*cos(theta); %y/(b/2);
a0 = a0_r - (a0_r - a0_t).*cos(theta); %y/(b/2);
c = c_r - (c_r - c_t).*cos(theta); %y/(b/2);
B = geo_AOA - AOA_0L; %B in Ax=B equation
for j = 1:N
    n(j) = (2*j - 1);
end
for i = 1:N
    for j = 1:N
        A(i,j) = (4*b/(a0(i)*c(i)))*sin(n(j)*theta(i)) + (n(j)*((sin(n(j)*theta(i)))/sin(theta(i)))); %A in Ax=b equation
    end
end
B = B';
An = A\B;

% x = linspace(0,b,numel(c));
% S = trapz(x,c);
S = ((c_r + c_t)/2)*b;
AR = b^2/S;
c_L = pi*AR*An(1);
delta = [];
for i = 2:N
    delta = [delta n(i)*(An(i)/An(1))^2];
end
delta = sum(delta);
e = 1/(1 + delta);
c_Di = c_L^2/(pi*e*AR);

%Calculate lift and induced drag for freestream velocity of 150 mph at sea
%level
% v_inf = 150;
% v_inf = convvel(150,'mph','ft/s');
% rho = .002378; %slug/ft^3
% q_inf = .5*rho*v_inf^2;
% L = c_L*S*q_inf;
% D_i = c_Di*S*q_inf;
% fprintf('The lift for the wing at a velocity of 150 miles/hour at sea level is %.3f pounds\n',L)
% fprintf('The induced drag for the wing at a velocity of 150 miles/hour at sea level is %.3f pounds\n',D_i)
%% Question 2 - find number of odd terms for requested relative error amounts, 5%,1%,.1%

end





