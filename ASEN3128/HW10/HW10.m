clear all;close all;clc
%define constants
global theta_0 u0 A
g = 9.81;
S = 510.9667;%m^2
b = 59.6433;%m
c_bar = 8.3241;%m
theta_0 = 0;
W = 2.83176e6;%Newtons
m = W/9.81;
V = 157.886;%m/s
rho = .6532;%kg/m^3
u0 = V;%m/s
xi = deg2rad(-6.8);
C_w0 = W/(.5*rho*(u0^2)*S);%weight coeff
I_x = 24675586.69;%kgm^2
I_y = 44877574.145;%kgm^2
I_z = 67384152.115;%kgm^2
I_zx = 1315143.4;%kgm^2
I_xs = I_x*cos(xi)^2 + I_z*sin(xi)^2 + I_zx*sin(2*xi);
I_zs = I_x*sin(xi)^2 + I_z*cos(xi)^2 - I_zx*sin(2*xi);
I_zxs = -1*(.5*(I_x - I_z)*sin(2*xi) + I_zx*(sin(xi)^2 - cos(xi)^2));
I_xp = (I_xs*I_zs - I_zxs^2)/I_zs;
I_zp = (I_xs*I_zs - I_zxs^2)/I_xs;
I_zxp = I_zxs/(I_xs*I_zs - I_zxs^2);
%non-dim stability matrix
NDvals = [-.8771 -.2797 .1946;...
              0 -.3295 -.04073;...
              0 .304 -.2737];
Lat_dim_derivs = [.5*rho*u0*S*NDvals(1,1) .5*rho*u0*b*S*NDvals(1,2)...
                    .5*rho*u0*b*S*NDvals(1,3);...
               .25*rho*u0*b*S*NDvals(2,1) .25*rho*u0*b^2*S*NDvals(2,2)...
               .25*rho*u0*b^2*S*NDvals(2,3);
               .25*rho*u0*b*S*NDvals(3,1) .25*rho*u0*b^2*S*NDvals(3,2)...
               .25*rho*u0*b^2*S*NDvals(3,3);];
%Conversion from body frame to stability frame
stab_lat(1,1) = Lat_dim_derivs(1,1);
stab_lat(2,1) = Lat_dim_derivs(2,1)*cos(xi) - Lat_dim_derivs(3,1)*sin(xi);
stab_lat(3,1) = Lat_dim_derivs(3,1)*cos(xi) + Lat_dim_derivs(2,1)*sin(xi);
stab_lat(1,2) = Lat_dim_derivs(1,2)*cos(xi) - Lat_dim_derivs(1,3)*sin(xi);
stab_lat(2,2) = Lat_dim_derivs(2,2)*cos(xi)^2 - (Lat_dim_derivs(3,2) + Lat_dim_derivs(2,3))*sin(xi)*cos(xi) + Lat_dim_derivs(3,3)*sin(xi)^2;
stab_lat(3,2) = Lat_dim_derivs(3,2)*cos(xi)^2 - (Lat_dim_derivs(3,3) - Lat_dim_derivs(2,2))*sin(xi)*cos(xi) - Lat_dim_derivs(2,3)*sin(xi)^2;
stab_lat(1,3) = Lat_dim_derivs(1,3)*cos(xi) + Lat_dim_derivs(1,2)*sin(xi);
stab_lat(2,3) = Lat_dim_derivs(2,3)*cos(xi)^2 - (Lat_dim_derivs(3,3) - Lat_dim_derivs(2,2))*sin(xi)*cos(xi) - Lat_dim_derivs(3,2)*sin(xi)^2;
stab_lat(3,3) = Lat_dim_derivs(3,3)*cos(xi)^2 + (Lat_dim_derivs(3,2) + Lat_dim_derivs(2,3))*sin(xi)*cos(xi) + Lat_dim_derivs(2,2)*sin(xi)^2;
%Construct the A matrix from the values above after having been converted
%to the stability frame
A(1,1) = stab_lat(1,1)/m;
A(2,1) = (stab_lat(1,2)/I_xp) + I_zxp*stab_lat(1,3);
A(3,1) = I_zxp*stab_lat(1,2) + (stab_lat(1,3)/I_zp);
A(4,1) = 0;
A(1,2) = stab_lat(2,1)/m;
A(2,2) = (stab_lat(2,2)/I_xp) + I_zxp*stab_lat(2,3);
A(3,2) = I_zxp*stab_lat(2,2) + (stab_lat(2,3)/I_zp);
A(4,2) = 1;
A(1,3) = (stab_lat(3,1)/m) - u0;
A(2,3) = (stab_lat(3,2)/I_xp) + I_zxp*stab_lat(3,3);
A(3,3) = I_zxp*stab_lat(3,2) + (stab_lat(3,3)/I_zp);
A(4,3) = tan(theta_0);
A(1,4) = g*cos(theta_0);
A(2,4) = 0;
A(3,4) = 0;
A(4,4) = 0;
%Determine eigen values of A matrix
[V,D] = eig(A);
%observing the A matrix eigen values in D, the first two (complex conj.),
%represents the Dutch roll mode, the larger negative real value is for the
%roll mode, and the smaller negative real value is the spiral mode. 
%determine time constants for each using TC = 1/(wn*zeta)
%Dutch Roll time constant
DRwn = sqrt(real(D(1,1))^2 + imag(D(1,1))^2);
DRzeta = -real(D(1,1))/DRwn;
DRTC = 1/(DRzeta*DRwn);
%Roll time constant
Rwn = sqrt(real(D(3,3))^2 + imag(D(3,3))^2);
Rzeta = -real(D(3,3))/Rwn;
RTC = 1/(Rwn*Rzeta);
%Spiral time constant
Swn = sqrt(real(D(4,4))^2 + imag(D(4,4))^2);
Szeta = -real(D(4,4))/Swn;
STC = 1/(Szeta*Swn);
%Question 4
%pull fancy yv,nr, and nv from A matrix
fYv = A(1,1);
fNv = A(3,1);
fNr = A(3,3);
%characteristic eqn. for Dutch roll eig. val. approx.
syms lambda
DRapprox = lambda^2 - (fYv + fNr)*lambda + (fYv*fNr +u0*fNv) == 0;
DReig = solve(DRapprox,lambda)
%Question 5a
tspan = [0:100];
inish_condish = [0 0 0 0 0 0];
[t,y] = ode45('hw10ode',tspan,inish_condish);
figure(1)
suptitle('Trim State')
subplot(3,2,1)
plot(t,y(:,1))
title('\Delta V')
xlabel('Time (s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta p')
xlabel('Time(s)')
ylabel('\Delta p (m/s)')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta r')
xlabel('Time(s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \phi')
xlabel('Time(s)')
ylabel('\Delta \phi (rads)')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta \psi')
xlabel('Time(s)')
ylabel('\Delta \psi (rads)')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Y')
xlabel('Time(s)')
ylabel('\Delta Y (m)')


%Question 5bi
tspan = [0:100];
inish_condish = [10 0 0 0 0 0];
[t,y] = ode45('hw10ode',tspan,inish_condish);
figure(2)
suptitle('\Delta v(0) = 10 m/s')
subplot(3,2,1)
plot(t,y(:,1))
title('\Delta V')
xlabel('Time(s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta P')
xlabel('Time(s)')
ylabel('\Delta p (m/s)')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta r')
xlabel('Time(s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \phi')
xlabel('Time(s)')
ylabel('\Delta \phi (rads)')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta \psi')
xlabel('Time(s)')
ylabel('\Delta \psi (rads)')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Y')
xlabel('Time(s)')
ylabel('\Delta Y (m)')
%Question 5bii
tspan = [0:100];
inish_condish = [0 0.15 0 0 0 0];
[t,y] = ode45('hw10ode',tspan,inish_condish);
figure(3)
suptitle('\Delta p(0) = 0.15 rads/s')
subplot(3,2,1)
plot(t,y(:,1))
title('\Delta V')
xlabel('Time(s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta P')
xlabel('Time(s)')
ylabel('\Delta p (m/s)')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta r')
xlabel('Time(s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \phi')
xlabel('Time(s)')
ylabel('\Delta \phi (rads)')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta \psi')
xlabel('Time(s)')
ylabel('\Delta \psi (rads)')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Y')
xlabel('Time(s)')
ylabel('\Delta Y (m)')
%Question 5biii
tspan = [0:100];
inish_condish = [-1.8563 -0.4185 0.0311 0.6148 0 0];
[t,y] = ode45('hw10ode',tspan,inish_condish);
figure(4)
suptitle('Case iii initial conditions')
subplot(3,2,1)
plot(t,y(:,1))
title('\Delta V')
xlabel('Time(s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta P')
xlabel('Time(s)')
ylabel('\Delta p (m/s)')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta r')
xlabel('Time(s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \phi')
xlabel('Time(s)')
ylabel('\Delta \phi (rads)')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta \psi')
xlabel('Time(s)')
ylabel('\Delta \psi (rads)')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Y')
xlabel('Time(s)')
ylabel('\Delta Y (m)')
%Question 5biv
tspan = [0:100];
inish_condish = [2.9477 -0.0063 0.0758 1.2431 0 0];
[t,y] = ode45('hw10ode',tspan,inish_condish);
figure(5)
suptitle('Case iv initial conditions')
subplot(3,2,1)
plot(t,y(:,1))
title('\Delta V')
xlabel('Time(s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
title('\Delta P')
xlabel('Time(s)')
ylabel('\Delta p (m/s)')
subplot(3,2,3)
plot(t,y(:,3))
title('\Delta r')
xlabel('Time(s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
title('\Delta \phi')
xlabel('Time(s)')
ylabel('\Delta \phi (rads)')
subplot(3,2,5)
plot(t,y(:,5))
title('\Delta \psi')
xlabel('Time(s)')
ylabel('\Delta \psi (rads)')
subplot(3,2,6)
plot(t,y(:,6))
title('\Delta Y')
xlabel('Time(s)')
ylabel('\Delta Y (m)')




