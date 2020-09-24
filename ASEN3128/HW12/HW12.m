clear all;close all;clc

%define constants
global theta_0 u0 A Alat_aug Blat_aug K
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
%ND control derivatives
NDcont = [0 -1.368e-2 -1.973e-4;0.1146 6.976e-3 -0.1257];
%Dimensional control derivatives
DMC(1,1) = NDcont(1,1)*.5*rho*u0^2*S;
DMC(1,2) = NDcont(1,2)*.5*rho*u0^2*S*b;
DMC(1,3) = NDcont(1,3)*.5*rho*u0^2*S*b;
DMC(2,1) = NDcont(2,1)*.5*rho*u0^2*S;
DMC(2,2) = NDcont(2,2)*.5*rho*u0^2*S*b;
DMC(2,3) = NDcont(2,3)*.5*rho*u0^2*S*b;
%construct 4X2 lateral B matrix
B(1,1) = DMC(1,1)/m;
B(1,2) = DMC(2,1)/m;
B(2,1) = (DMC(1,2)/I_xp) + I_zxp*DMC(1,3);
B(2,2) = (DMC(2,2)/I_xp) + I_zxp*DMC(2,3);
B(3,1) = I_zxp*DMC(1,2) + (DMC(1,3)/I_zp);
B(3,2) = I_zxp*DMC(2,2) + (DMC(2,3)/I_zp);
B(4,1) = 0;
B(4,2) = 0;
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
%determine E and D varibles to compute real root for spiral mode, -E/D
%define fancy letters used to calculate E and D
fyv = A(1,1); %fancy y_v
flv = A(2,1); %fancy l_v
fnv = A(3,1); %fancy n_v
flp = A(2,2); %fancy l_p
fnp = A(3,2); %fancy n_p
fyr = A(1,3); %fancy y_r
flr = A(2,3); %fancy l_r
fnr = A(3,3); %fancy n_r
E = g*((flv*fnr - flr*fnv)*cos(theta_0) + (flp*fnv - flv*fnp)*sin(theta_0));
D = -g*(flv*cos(theta_0) + fnv*sin(theta_0)) + u0*(flv*fnp - flp*fnv);
%calculate spiral mode eigenvalue approximation
spiral_lambda = -E/D;
%determine roll mode eigenvalue approx
roll_lambda = flp;
Alat_aug = zeros(6,6);
Alat_aug((1:4),(1:4)) = A;
Alat_aug(5,(1:6)) = [0 0 sec(theta_0) 0 0 0];   
Alat_aug(6,:) = [1 0 0 0 u0*cos(theta_0) 0];
Blat_aug = zeros(6,2);
Blat_aug((1:4),:) = B;
Blat_aug((5:6),:) = 0;
K = [0 1 0 0 0 0;0 0 1.5 0 0 0];
tspan = [0:100];
inish_condish = [10 -.14 0.05 0 0 0];
[t,y] = ode45('hw12ode',tspan,inish_condish);
y(:,4) = y(:,4)*(180/pi);
y(:,5) = y(:,5)*(180/pi);
% %% Question 3
% %3a
k1 = [0:0.5:5];
k2 = [0:0.5:5];
K = [0 k1(1) 0 0 0 0;0 0 k2(1) 0 0 0];
%Kopt = [0 1 0 0 0 0;0 0 1.5 0 0 0]
figure
hold on
grid on
grid minor
Alat_aug_cl = Alat_aug + Blat_aug*K;
[V,D] = eig(Alat_aug_cl);
plot(diag(D),'g*','MarkerSize',15)
K = [0 k1(end) 0 0 0 0; 0 0 k2(end) 0 0 0];
Alat_aug_cl = Alat_aug + Blat_aug*K;
[V,D] = eig(Alat_aug_cl);
plot(diag(D),'r*','MarkerSize',15)
% legend('First','Last')
counter = 1;
for i = 2:length(k1) - 1
    for j = 2:length(k2) - 1
        K = [0 k1(i) 0 0 0 0;0 0 k2(i) 0 0 0];
        Alat_aug_cl = Alat_aug + Blat_aug*K;
        [V,D] = eig(Alat_aug_cl);
        plot(diag(D),'b o')
        spiral_tc(counter) = -1/real(D(5,5));
        roll_tc(counter) = -1/real(D(6,6));
        DR_wn = norm(D(3,3));
        DR_zeta = -real(D(3,3))/DR_wn;
        DR_tc(counter) = 1/(DR_zeta*DR_wn);
        counter = counter + 1;
    end
end
title('Eigenvalues for varying K_{ap} and K_{rr}','FontSize',15)
xlabel('Real','FontSize',15)
ylabel('Imgaginary','FontSize',15)
legend('Min K value','Max K value')
figure(2)
hold on
grid on 
grid minor
plot(spiral_tc)
plot(roll_tc)
plot(DR_tc)
title('Time Constant Vs K','FontSize',15)
xlabel('K','FontSize',18)
ylabel('Time Constant (s)','FontSize',15)
legend('Spiral Time Constant','Roll Time Constant','Dutch Roll Time Constant')
K = [0 1 0 0 0 0;0 0 1.5 0 0 0];
Alat_aug_cl = Alat_aug + Blat_aug*K;
% saveas(gcf,'3a.jpg')
for i = 1:length(y)
    deflect(:,i) = K*y(i,:)';
end
deflect = deflect';
deflect = deflect*(180/pi);
figure
hold on
grid on
grid minor
plot(t,deflect)
title('Aileron and Rudder Deflection Vs Time','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Degrees','FontSize',15)
legend('\Delta\delta a','\Delta\delta r')
[V,D] = eig(Alat_aug_cl);
spiral_tc = -1/real(D(5,5));
roll_tc = -1/real(D(6,6));
DR_wn = norm(D(3,3));
DR_zeta = -real(D(3,3))/DR_wn;
DR_tc = 1/(DR_zeta*DR_wn);
%feed back state variables 
figure
suptitle('State Variables Vs Time With Feedback')
subplot(3,2,1)
plot(t,y(:,1))
xlabel('Time (s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
xlabel('Time (s)')
ylabel('\Delta p (rads/s)')
subplot(3,2,3)
plot(t,y(:,3))
xlabel('Time (s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
xlabel('Time (s)')
ylabel('\Delta \phi (degrees)')
subplot(3,2,5)
plot(t,y(:,5))
xlabel('Time (s)')
ylabel('\Delta \psi (degrees)')
subplot(3,2,6)
plot(t,y(:,6))
xlabel('Time (s)')
ylabel('Meters')
%open loop state variable response
K = [0 0 0 0 0 0;0 0 0 0 0 0];
tspan = [0:100];
inish_condish = [10 -.14 0.05 0 0 0];
[t,y] = ode45('hw12ode',tspan,inish_condish);
y(:,4) = y(:,4)*(180/pi);
y(:,5) = y(:,5)*(180/pi);
figure
suptitle('State Variables Vs Time Without Feedback')
subplot(3,2,1)
plot(t,y(:,1))
xlabel('Time (s)')
ylabel('\Delta V (m/s)')
subplot(3,2,2)
plot(t,y(:,2))
xlabel('Time (s)')
ylabel('\Delta p (rads/s)')
subplot(3,2,3)
plot(t,y(:,3))
xlabel('Time (s)')
ylabel('\Delta r (rads/s)')
subplot(3,2,4)
plot(t,y(:,4))
xlabel('Time (s)')
ylabel('\Delta \phi (degrees)')
subplot(3,2,5)
plot(t,y(:,5))
xlabel('Time (s)')
ylabel('\Delta \psi (degrees)')
subplot(3,2,6)
plot(t,y(:,6))
xlabel('Time (s)')
ylabel('Meters')

