%{
Midterm Project -- Reconstruction of Galileo's interplanetary cruise to
Jupiter. Construct PCPs using given baseline mission dates and then create
full mission scenario in GMAT for analysis and visualization. 
%}
clear all;close all;clc
%% Part 1
% Construct PCPs for each segment ofthe mission and pick point in each PCP
% that corresponds to Galileo's traj.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Launch-Venus PCP

%Defien specified event dates
Event.LaunchJD = 2447807.5;
Event.VGAJD =  2447932.5;
Event.EGA1JD = 2448235.5;
Event.EGA2JD = 2448965.984378;
Event.JOIJD = 2450164;

mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.LaunchJD - 20:Event.LaunchJD + 70];
Arr_Dates = [Event.VGAJD - 20:Event.VGAJD + 70];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Venus,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Venus(i,j) = norm(V_final - V_Venus);
        C3(i,j) = (norm(V_initial - V_Earth))^2;
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [7.3 7.5 7.8 8 9 13 18 19 20 21 25];
V_inf_contours = [1 2 2.5 3 3.2 3.8 4 4.5 5 6 7 8 12];
dt_contours = [50 55 60 65 70 75 80 90 100 110 120 130 140 150 160 170 180 190];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,C3,C3_contours,'r','LineWidth',1);
clabel(cs1,h1,'FontSize',16)
[cs2,h2] = contour(x_vals,y_vals,vinf_Venus,V_inf_contours,'b','LineWidth',1);
clabel(cs2,h2,'FontSize',16)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k','LineWidth',1);
clabel(cs3,h3,'FontSize',16)
plot(20,20,'g*','MarkerSize',15)
set(gca,'FontSize',20)
xlabel('Days Past 18 Sept 1989')
ylabel('Days Past 21 Jan 1990')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Venus}$,$\frac{km}{s}$','TOF, days','Galileo Traj.','FontSize',18,'Interpreter','latex')
title('Earth-Venus Galileo Pork Chop Plot','FontSize',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Venus-EGA1 PCP

mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.VGAJD - 120:Event.VGAJD + 20];
Arr_Dates = [Event.EGA1JD - 250:Event.EGA1JD + 20];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Venus,V_Venus] = calcposvel(EphemD.Venus.a,EphemD.Venus.e,EphemD.Venus.i,EphemD.Venus.Omega,EphemD.Venus.w,EphemD.Venus.nu,mu);
        [R_Earth,V_Earth] = calcposvel(EphemA.Earth.a,EphemA.Earth.e,EphemA.Earth.i,EphemA.Earth.Omega,EphemA.Earth.w,EphemA.Earth.nu,mu);
        [V_initial,V_final] = solvelambert(R_Venus,R_Earth,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Venus(i,j) = norm(V_initial - V_Venus);
        vinf_Earth(i,j) = norm(V_final - V_Earth);
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

V_inf_Venus_contours = [2.5 2.8 3 3.5 4 4.5 5 6 7 11 14];
V_inf_Earth_contours = [2.5 2.8 3 3.5 4 4.5 5 6 6.5 7 8 10 13 15];
dt_contours = [200 220 240 260 280 300 320 340 360 380 400 420 440];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,vinf_Venus,V_inf_Venus_contours,'r','LineWidth',1);
clabel(cs1,h1,'FontSize',16)
[cs2,h2] = contour(x_vals,y_vals,vinf_Earth,V_inf_Earth_contours,'b','LineWidth',1);
clabel(cs2,h2,'FontSize',16)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k','LineWidth',1);
clabel(cs3,h3,'FontSize',16)
plot(120,250,'g*','MarkerSize',15)
set(gca,'FontSize',20)
xlabel('Days Past 13 Oct 1989')
ylabel('Days Past 4 Apr 1990')
legend('$V_{\infty,Venus}$,$\frac{km}{s}$','$V_{\infty,Earth}$,$\frac{km}{s}$','TOF, days','Galileo Traj.','FontSize',18,'Interpreter','latex')
title('Venus-EGA1 Galileo Pork Chop Plot','FontSize',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EGA2-JOI PCP
clear vinf_Earth dt
mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.EGA2JD - 120:Event.EGA2JD + 120];
Arr_Dates = [Event.JOIJD - 300:Event.JOIJD + 300];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Jupiter,V_Jupiter] = calcposvel(EphemA.Jupiter.a,EphemA.Jupiter.e,EphemA.Jupiter.i,EphemA.Jupiter.Omega,EphemA.Jupiter.w,EphemA.Jupiter.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Jupiter,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Earth(i,j) = norm(V_initial - V_Earth);
        vinf_Jupiter(i,j) = norm(V_final - V_Jupiter);
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

V_inf_Earth_contours = [8.5 8.8 9 9.5 10 11 13 15];
V_inf_Jupiter_contours = [5 5.3 5.5 5.8 6 6.3 6.8 7 7.5 8 9 10 12 15];
dt_contours = [750 780 800 850 900 950 1001 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1660 1650];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,vinf_Earth,V_inf_Earth_contours,'r','LineWidth',1);
clabel(cs1,h1,'FontSize',16)
[cs2,h2] = contour(x_vals,y_vals,vinf_Jupiter,V_inf_Jupiter_contours,'b','LineWidth',1);
clabel(cs2,h2,'FontSize',16)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k','LineWidth',1);
clabel(cs3,h3,'FontSize',16)
plot(120,300,'g*','MarkerSize',15)
set(gca,'FontSize',20)
xlabel('Days Past 13 Oct 1989')
ylabel('Days Past 4 Apr 1990')
legend('$V_{\infty,Earth}$,$\frac{km}{s}$','$V_{\infty,Jupiter}$,$\frac{km}{s}$','TOF, days','Galileo Traj.','FontSize',18,'Interpreter','latex')
title('EGA2-JOI Galileo Pork Chop Plot','FontSize',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2 
% Construct Earth-Earth resonant orbit plot and state phi value that was
% selected and why

mu_S = 1.32712440018e11;
mu_E = 3.98600433e5;
mu_V = 3.24858599e5;
AU = 1.49597870700e8;
mu_J = 1.266865361e8;

Re = 6378.14;
Rv = 6051.8;
Rj = 71492;

% Event.Launch = 2458951.5;
% Event.VGA1 = 2459121.5;
% Event.EGA1 = 2459438.5;
% Event.EGA2 = 2460534.226567;
% Event.JOI = 2461300.75;

% Calculate hyperbolic excess velocity of S/C upon arrivalt at EGA1
VGA1Ephem = Meeus(Event.VGAJD);
EGA1Ephem = Meeus(Event.EGA1JD);

[R_VGA1,Vp_VGA1] = calcposvel(VGA1Ephem.Venus.a,VGA1Ephem.Venus.e,VGA1Ephem.Venus.i,VGA1Ephem.Venus.Omega,VGA1Ephem.Venus.w,VGA1Ephem.Venus.nu,mu_S);
[R_EGA1,Vp_EGA1] = calcposvel(EGA1Ephem.Earth.a,EGA1Ephem.Earth.e,EGA1Ephem.Earth.i,EGA1Ephem.Earth.Omega,EGA1Ephem.Earth.w,EGA1Ephem.Earth.nu,mu_S);
% V_inf_in = vf - v0V, V_inf_out = v0_V_FB - v0V;
[v_sc_VGA1,v_sc_EGA1] = solvelambert(R_VGA1,R_EGA1,(Event.EGA1JD - Event.VGAJD)*86400,0);

V_inf_in_EGA1 = v_sc_EGA1 - Vp_EGA1;

%Calculate hyperbolic excess velocity of S/C up departure at EGA1
EGA2Ephem = Meeus(Event.EGA2JD);
JOIEphem = Meeus(Event.JOIJD);

[R_EGA2,Vp_EGA2] = calcposvel(EGA2Ephem.Earth.a,EGA2Ephem.Earth.e,EGA2Ephem.Earth.i,EGA2Ephem.Earth.Omega,EGA2Ephem.Earth.w,EGA2Ephem.Earth.nu,mu_S);
[R_JOI,Vp_JOI] = calcposvel(JOIEphem.Jupiter.a,JOIEphem.Jupiter.e,JOIEphem.Jupiter.i,JOIEphem.Jupiter.Omega,JOIEphem.Jupiter.w,JOIEphem.Jupiter.nu,mu_S);

[v_sc_EGA2,v_sc_JOI] = solvelambert(R_EGA2,R_JOI,(Event.JOIJD - Event.EGA2JD)*86400,0);

V_inf_out_EGA2 = v_sc_EGA2 - Vp_EGA2;

V_inf_in_JOI = v_sc_JOI - Vp_JOI;

%Construct periapse radius of each gravity assist as func. of phi
%first flyby
P_Earth = 2*365.242189*86400;
a = ((P_Earth/(2*pi))^2*mu_S)^(1/3);
V_sc_post_EGA1 = sqrt(mu_S*((2/norm(R_EGA1)) - (1/a)));
V_inf = (norm(V_inf_out_EGA2) + norm(V_inf_in_EGA1))/2;
V_sc_sun = sqrt((2*mu_S)/norm(R_EGA1) - (mu_S/a));

num = -norm(V_sc_sun)^2 + V_inf^2 + norm(Vp_EGA1)^2;
denom = 2*V_inf*norm(Vp_EGA1);
theta = acos((num/denom));
%construct VNC frame transformation
V_hat = Vp_EGA1/norm(Vp_EGA1);
N_hat = (cross(R_EGA1,Vp_EGA1))/norm(cross(R_EGA1,Vp_EGA1));
C_hat = cross(V_hat,N_hat);
T = [V_hat N_hat C_hat];

phi = [0:0.01:2*pi];

for i = 1:numel(phi)
    V_inf_out_EGA1 = V_inf*[cos(pi - theta);sin(pi - theta)*cos(phi(i));-sin(pi - theta)*sin(phi(i))];
    V_inf_out_EGA1(:,i) = T*V_inf_out_EGA1;
    V_inf_in_EGA2(:,i) = V_inf_out_EGA1(:,i) + Vp_EGA1 - Vp_EGA2;
    Psi_EGA1(i) = acos((dot(V_inf_in_EGA1,V_inf_out_EGA1(:,i)))/(norm(V_inf_in_EGA1)*norm(V_inf_out_EGA1(:,i))));
    Psi_EGA2(i) = acos((dot(V_inf_in_EGA2(:,i),V_inf_out_EGA2))/(norm(V_inf_in_EGA2(:,i))*norm(V_inf_out_EGA2)));
    rp_EGA1(i) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA1(i))/2) - 1));
    rp_EGA2(i) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA2(i))/2) - 1));
end
phiEGA1 = find(rp_EGA1 > (Re+300));
phiEGA2 = find(rp_EGA2 > (Re+300));
acceptable_phi = intersect(phiEGA1,phiEGA2);
plot_acc_phix = [rad2deg(phi(min(acceptable_phi))) rad2deg(phi(min(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) ];
plot_acc_phiy = [0 12000 12000 0];    

rplim = (Re+300)*ones(numel(phi),1);

figure
hold on
plot(phi*(180/pi),rp_EGA1,'LineWidth',1)
plot(phi*(180/pi),rp_EGA2,'LineWidth',1)
plot(phi*(180/pi),rplim,'--k','LineWidth',1)
patch(plot_acc_phix,plot_acc_phiy,'gree','FaceAlpha',0.2)
set(gca,'FontSize',20)
xlim([0 360])
ylim([0 12000])
grid on
grid minor

xlabel('\phi (Degrees)')
ylabel('r_p (km)')
title('Periapse Radius r_p VS \phi','FontSize',24)
legend('EGA1 r_p','EGA2 r_p','Minimum r_p','Acceptable \phi region','FontSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3
% Use code to fill in theoretical data for the Galileo traj. in table 2

%Launch Params -- C3, RLA, DLA
EphemD = Meeus(Event.LaunchJD);
EphemA = Meeus(Event.VGAJD);
[R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
[R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
[V_initial,V_final] = solvelambert(R_Earth,R_Venus,(Event.VGAJD - Event.LaunchJD)*86400,0);
vinf_Venus = norm(V_final - V_Venus);
C3_Launch = (norm(V_initial - V_Earth))^2;

V_inf_out_launch = V_initial - V_Earth;
RLA = rad2deg(atan2(V_inf_out_launch(2),V_inf_out_launch(1)));
DLA = rad2deg(asin(V_inf_out_launch(3)/norm(V_inf_out_launch)));
%radius and altitude of closest approach for each flyby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VGA
EphemD = Meeus(Event.LaunchJD);
EphemA = Meeus(Event.VGAJD);
[R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
[R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
[V_initial,V_final] = solvelambert(R_Earth,R_Venus,(Event.VGAJD - Event.LaunchJD)*86400,0);
vinf_Venus_in = V_final - V_Venus;
vinf_Earth = V_initial - V_Earth;

EphemD = Meeus(Event.VGAJD);
EphemA = Meeus(Event.EGA1JD);
[R_Venus,V_Venus] = calcposvel(EphemD.Venus.a,EphemD.Venus.e,EphemD.Venus.i,EphemD.Venus.Omega,EphemD.Venus.w,EphemD.Venus.nu,mu);
[R_Earth,V_Earth] = calcposvel(EphemA.Earth.a,EphemA.Earth.e,EphemA.Earth.i,EphemA.Earth.Omega,EphemA.Earth.w,EphemA.Earth.nu,mu);

R_EGA1 = R_Earth;
Vp_EGA1 = V_Earth;

[V_initial,V_final] = solvelambert(R_Venus,R_Earth,(Event.EGA1JD - Event.VGAJD)*86400,0);
vinf_Venus_out = V_initial - V_Venus;
vinf_Earth_in = V_final - V_Earth;
V_inf_in_EGA1 = vinf_Earth_in;

V_inf_in = vinf_Venus_in;
V_inf_out = vinf_Venus_out;

V_inf = norm(V_inf_in);
% mu_J = 1.266865361e8; %km^3/s^2
% R_Jupiter = 71492; %Jupiter radius, km

k_hat = [0 0 1];
S_hat = V_inf_in/norm(V_inf_in);
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat);
h_hat = cross(V_inf_in,V_inf_out)/norm(cross(V_inf_in,V_inf_out));
B_hat = cross(S_hat,h_hat);

a = -mu_V/V_inf^2;
num = dot(V_inf_out,V_inf_in);
denom = norm(V_inf_out)*norm(V_inf_in);
Psi = acos(num/denom);
rho = (pi - Psi)/2;
rp = (mu_V/V_inf^2)*((1/cos(rho))-1);
c = a/cos(rho);
b = sqrt(a^2*((c/a)^2 - 1));
B = b*B_hat;

theta = acos((dot(B_hat,T_hat))/(norm(B_hat)*norm(T_hat)));
B_T = dot(B_hat,T_hat);
B_R = dot(B_hat,R_hat);
normB = norm(B);

num = dot(V_inf_in,V_inf_out);
denom = norm(V_inf_in)*norm(V_inf_out);
Psi = acos(num/denom);

rp = (mu_V/norm(V_inf_in)^2)*((1/(cos((pi-Psi)/2))) - 1);

flyby_alt = rp - Rv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EGA1/EGA2
%compute EGA2 to JOI params in order to get EGA1 to EGA2 params/ vinfs
EGA2Ephem = Meeus(Event.EGA2JD);
JOIEphem = Meeus(Event.JOIJD);

[R_EGA2,Vp_EGA2] = calcposvel(EGA2Ephem.Earth.a,EGA2Ephem.Earth.e,EGA2Ephem.Earth.i,EGA2Ephem.Earth.Omega,EGA2Ephem.Earth.w,EGA2Ephem.Earth.nu,mu);
[R_JOI,Vp_JOI] = calcposvel(JOIEphem.Jupiter.a,JOIEphem.Jupiter.e,JOIEphem.Jupiter.i,JOIEphem.Jupiter.Omega,JOIEphem.Jupiter.w,JOIEphem.Jupiter.nu,mu);

[v_sc_EGA2,v_sc_JOI] = solvelambert(R_EGA2,R_JOI,(Event.JOIJD - Event.EGA2JD)*86400,0);

V_inf_out_EGA2 = v_sc_EGA2 - Vp_EGA2;

%Construct periapse radius of gravity assist as func. of phi
%first flyby
P_Earth = 2*365.242189*86400;
a = ((P_Earth/(2*pi))^2*mu_S)^(1/3);
V_sc_post_EGA1 = sqrt(mu_S*((2/norm(R_EGA1)) - (1/a)));
% V_inf = (norm(V_inf_out_EGA1) + norm(V_inf_in_EGA1))/2;
V_inf = norm(V_inf_in_EGA1);
V_sc_sun = sqrt((2*mu_S)/norm(R_EGA1) - (mu_S/a));
phi = deg2rad(100);
% phi = deg2rad(99.0990990990991);

num = -norm(V_sc_sun)^2 + V_inf^2 + norm(Vp_EGA1)^2;
denom = 2*V_inf*norm(Vp_EGA1);
theta = acos((num/denom));

V_hat = Vp_EGA1/norm(Vp_EGA1);
N_hat = (cross(R_EGA1,Vp_EGA1))/norm(cross(R_EGA1,Vp_EGA1));
C_hat = cross(V_hat,N_hat);
T = [V_hat N_hat C_hat];

V_inf_out_EGA1 = V_inf*[cos(pi - theta);sin(pi - theta)*cos(phi);-sin(pi - theta)*sin(phi)];
V_inf_out_EGA1 = T*V_inf_out_EGA1;
V_inf_in_EGA2 = V_inf_out_EGA1 + Vp_EGA1 - Vp_EGA2;
Psi_EGA1 = acos((dot(V_inf_in_EGA1,V_inf_out_EGA1))/(norm(V_inf_in_EGA1)*norm(V_inf_out_EGA1)));
Psi_EGA2 = acos((dot(V_inf_in_EGA2,V_inf_out_EGA2))/(norm(V_inf_in_EGA2)*norm(V_inf_out_EGA2)));
rp_EGA1 = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA1)/2) - 1));
rp_EGA2 = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA2)/2) - 1));

V_inf = (norm(V_inf_in_EGA1) + norm(V_inf_out_EGA1))/2;

%B-plane parameter calculation -- EGA1 
k_hat = [0 0 1];
S_hat = V_inf_in_EGA1/norm(V_inf_in_EGA1);
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat);
h_hat = cross(V_inf_in_EGA1,V_inf_out_EGA1)/norm(cross(V_inf_in_EGA1,V_inf_out_EGA1));
B_hat = cross(S_hat,h_hat);

a = -mu_E/V_inf^2;
num = dot(V_inf_out_EGA1,V_inf_in_EGA1); 
denom = norm(V_inf_out_EGA1)*norm(V_inf_in_EGA1);
Psi = acos(num/denom);
rho = (pi - Psi)/2;
rp = (mu_E/V_inf^2)*((1/cos(rho))-1);
c = a/cos(rho);
b = sqrt(a^2*((c/a)^2 - 1));
B = b*B_hat;

theta = acos((dot(B_hat,T_hat))/(norm(B_hat)*norm(T_hat)));
B_T = dot(B_hat,T_hat);
B_R = dot(B_hat,R_hat);
normB = norm(B);

num = dot(V_inf_in_EGA1,V_inf_out_EGA1);
denom = norm(V_inf_in_EGA1)*norm(V_inf_out_EGA1);
Psi = acos(num/denom);

rp = (mu_E/norm(V_inf_in_EGA1)^2)*((1/(cos((pi-Psi)/2))) - 1);

flyby_alt = rp - Re;

%second flyby

V_inf = (norm(V_inf_in_EGA2) + norm(V_inf_out_EGA2))/2;

%B-plane parameter calculation -- EGA2
k_hat = [0 0 1];
S_hat = V_inf_in_EGA2/norm(V_inf_in_EGA2);
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat);
h_hat = cross(V_inf_in_EGA2,V_inf_out_EGA2)/norm(cross(V_inf_in_EGA2,V_inf_out_EGA2));
B_hat = cross(S_hat,h_hat);

a = -mu_E/V_inf^2;
num = dot(V_inf_out_EGA2,V_inf_in_EGA2); 
denom = norm(V_inf_out_EGA2)*norm(V_inf_in_EGA2);
Psi = acos(num/denom);
rho = (pi - Psi)/2;
rp = (mu_E/V_inf^2)*((1/cos(rho))-1);
c = a/cos(rho);
b = sqrt(a^2*((c/a)^2 - 1));
B = b*B_hat;

theta = acos((dot(B_hat,T_hat))/(norm(B_hat)*norm(T_hat)));
B_T = dot(B_hat,T_hat);
B_R = dot(B_hat,R_hat);
normB = norm(B);

num = dot(V_inf_in_EGA2,V_inf_out_EGA2);
denom = norm(V_inf_in_EGA2)*norm(V_inf_out_EGA2);
Psi = acos(num/denom);

rp = (mu_E/norm(V_inf_in_EGA2)^2)*((1/(cos((pi-Psi)/2))) - 1);

flyby_alt = rp - Re;
%data found from websites/papers
a_j = 9814688.8; %km
e_j = 0.971233;
inc_j = 5.3; %degrees

rp_j = a_j*(1 - e_j);

V_hyp = sqrt(((2*mu_J)/rp_j) + norm(V_inf_in_JOI)^2);
V_el = sqrt(((2*mu_J)/rp_j) - (mu_J/a_j));

dv = V_el - V_hyp;