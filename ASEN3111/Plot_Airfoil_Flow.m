%Author:Sage Herrin
%Created: 9/21/18
%SID:106071909
%Comp Lab 2
%ASEN 3111, Aerodynamics
function [V P] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
%%Write matlab function that plots stream lines, equipotential lines, and
%%pressure contours for flow about a thin symmetric airfoil
%Domain
xmin = -3*c/c;
xmax = 3*c/c;
ymin = -c;
ymax = c;
%Number of grid points
dx = 100;
dy = 100;
%Mesh over specified number of grid points
[x,y] = meshgrid(linspace(xmin,xmax,dx),linspace(ymin,ymax,dy));
delta_x = c/N;
airfoil_x = [delta_x/4:delta_x:c-delta_x/4];
airfoil_y = zeros(1,length(airfoil_x));
% num = 1 - (x/c);
% denom = x/c;
Psi = 0;
Phi = 0;
upper_g = 0;
%Convert angle to rads
alpha = deg2rad(alpha);
for i = 1:N
    lower_g = calc_g(airfoil_x(i),c,alpha,V_inf);
    upper_g = lower_g*delta_x;
    Psi = Psi + (upper_g/(2*pi))*log(radius(x,y,airfoil_x(i),airfoil_y(i)));
    Phi = Phi - ((upper_g/(2*pi)))*(atan2(-(y - airfoil_y(i)),-(x - airfoil_x(i))));
end
%Uniform flow components
u = V_inf*cos(alpha);
v = V_inf*sin(alpha);
%Combine uniform and vortex flow with superposition
Psi = Psi + u*y - v*x;
Phi = Phi + u*x + v*y;
%%Contours
%Psi Contour
% levmin = Psi(1,dx);
% levmax = Psi(dy,dx/2);
Psilevmin = min(min(Psi));
Psilevmax = max(max(Psi));
Psilevels = linspace(Psilevmin,Psilevmax,50);
%Phi Contour
Philevmin = min(min(Phi));
Philevmax = max(max(Phi));
Philevels = linspace(Philevmin,Philevmax,50);
%Countour plots
figure(1)
hold on
contourf(x,y,Psi,Psilevels)
h = plot([0 c],[0 0]);
set(h(1),'linewidth',2);
xlabel('X Position (Meters)')
ylabel('Y Position (Meters)')
title('StreamLines')
hold off
figure(2)
hold on
contourf(x,y,Phi,Philevels)
h = plot([0 c],[0 0]);
set(h(1),'linewidth',2);
xlabel('X Position (Meters)')
ylabel('Y Position (Meters)')
title('Equipotential Lines')
hold off
%Calculate pressure for pressure contour 
%Velocity in field using the gradient of Phi
V = gradient(Phi,(xmax - xmin)/dx,(ymax - ymin)/dy);
q = 5*rho_inf*(V_inf)^2;
P = q.*(1 - (V./V_inf).^2) + p_inf;
%Pressure contour
Plevmin = min(min(P));
Plevmax = max(max(P));
Plevels = linspace(Plevmin,Plevmax,50);
figure(3)
hold on
contourf(x,y,P,Plevels)
colorbar;
h = plot([0 c],[0 0]);
set(h(1),'linewidth',2);
xlabel('X Position (Meters)')
ylabel('Y Position (Meters)')
title('Pressure Contour')
hold off




