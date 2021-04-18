clear all;close all;clc
%% Part I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1 -- use Lambert code to determine launch parameters (C3,RLA,DLA)
%and JGA arrival conditions (V_inf,JGA)
mu_S = 1.32712440018e11;
Launch_JD = 2453755.29167; %Jan 9, 2006 19:00:00
JGA_JD = 2454159.73681; %Feb 28, 2007 05:41:00
PCE_JD = 2457217.99931; %July 14, 2015 11:59:00

% Calculate Ephemeris at Launch
LaunchEphem = Meeus(Launch_JD);
JGAEphem = Meeus(JGA_JD);
[R_Earth,V_Earth] = calcposvel(LaunchEphem.Earth.a,LaunchEphem.Earth.e,LaunchEphem.Earth.i,LaunchEphem.Earth.Omega,LaunchEphem.Earth.w,LaunchEphem.Earth.nu,mu_S);
[R_Jupiter,V_Jupiter] = calcposvel(JGAEphem.Jupiter.a,JGAEphem.Jupiter.e,JGAEphem.Jupiter.i,JGAEphem.Jupiter.Omega,JGAEphem.Jupiter.w,JGAEphem.Jupiter.nu,mu_S);

[V_sc_Launch,V_sc_JGA] = solvelambert(R_Earth,R_Jupiter,(JGA_JD - Launch_JD)*86400,0);

C3 = (norm(V_sc_Launch - V_Earth))^2;

V_inf_out_launch = V_sc_Launch - V_Earth;

RLA = atan2(V_inf_out_launch(2),V_inf_out_launch(1));
DLA = asin(V_inf_out_launch(3)/norm(V_inf_out_launch));

V_inf_in_JGA = V_sc_JGA - V_Jupiter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 2 -- use Lambert code to determine JGA dep. conditions
%(v_inf_out_JGA) and Pluto/Charon arrival conditions (V_inf_in_Pluto)

PCEEphem = Meeus(PCE_JD);
[R_Pluto,V_Pluto] = calcposvel(PCEEphem.Pluto.a,PCEEphem.Pluto.e,PCEEphem.Pluto.i,PCEEphem.Pluto.Omega,PCEEphem.Pluto.w,PCEEphem.Pluto.nu,mu_S);

[V_sc_JGA_out,V_sc_Pluto] = solvelambert(R_Jupiter,R_Pluto,(PCE_JD - JGA_JD)*86400,0);

V_inf_out_JGA = V_sc_JGA_out - V_Jupiter;
V_inf_in_PCE = V_sc_Pluto - V_Pluto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 3 -- Explain why values of |V_inf_in_JGA| and |V_inf_out_JGA| are
%same magnitude, what might cause differences (done on answer sheet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 4 -- Use gravity-flyby code to determine B_plane params, turn
%angle, and altitude of closest approach w/ Jupiter for JGA.
V_inf_in = V_inf_in_JGA;
V_inf_out = V_inf_out_JGA;

V_inf = norm(V_inf_in);
mu_J = 1.266865361e8; %km^3/s^2
R_Jupiter = 71492; %Jupiter radius, km

k_hat = [0 0 1];
S_hat = V_inf_in/norm(V_inf_in);
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat);
h_hat = cross(V_inf_in,V_inf_out)/norm(cross(V_inf_in,V_inf_out));
B_hat = cross(S_hat,h_hat);

a = -mu_J/V_inf^2;
num = dot(V_inf_out,V_inf_in);
denom = norm(V_inf_out)*norm(V_inf_in);
Psi = acos(num/denom);
rho = (pi - Psi)/2;
rp = (mu_J/V_inf^2)*((1/cos(rho))-1);
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

rp = (mu_J/norm(V_inf_in)^2)*((1/(cos((pi-Psi)/2))) - 1);

flyby_alt = rp - R_Jupiter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 5 -- Calculate Heliocentric dv provided by Jupiter gravity assist
Helio_dv = V_sc_JGA_out - V_sc_JGA;
%% Part II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCPs -- Construct 2 PCPs, one of available trajectories between launch
% and JGA, and one between JGA and PCE
%Part 1 -- Launch-JGA trjaectory possibilities
mu = 1.32712440018e11; %mu of sun
Dep_Dates = [2453714.5:2453794.5];
Arr_Dates = [2454129.5:2454239.5];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Jupiter,V_Jupiter] = calcposvel(EphemA.Jupiter.a,EphemA.Jupiter.e,EphemA.Jupiter.i,EphemA.Jupiter.Omega,EphemA.Jupiter.w,EphemA.Jupiter.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Jupiter,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Jupiter(i,j) = norm(V_final - V_Jupiter);
        C3(i,j) = (norm(V_initial - V_Earth))^2;
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [70,80,90,100,110,115,120,130,140,150,160,170,180,200,210,230,250,280,310,320,330,340,350];
V_inf_contours = [10,11,12,13,13.3,13.5,13.8,14.2,14.5,15,16,18,19,19.5,20,20.5,21,21.3,21.5,21.6,22,22.5,23];
dt_contours = [260,280,300,310,335,350,370,400,420,450,480,510,520,530,550,580];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,C3,C3_contours,'r');
clabel(cs1,h1)
[cs2,h2] = contour(x_vals,y_vals,vinf_Jupiter,V_inf_contours,'b');
clabel(cs2,h2)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k');
clabel(cs3,h3)
plot(40,30,'g*','MarkerSize',15)
xlabel('Days Past 10 Dec 2005')
ylabel('Days Past 29 Jan 2007')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Jupiter}$,$\frac{km}{s}$','TOF, days','NH Traj. Dates','Interpreter','latex')
title('Earth-Jupiter New Horizons Pork Chop Plot')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 2 -- JGA-PCE PCP for possible trajectories
mu = 1.32712440018e11; %mu of sun
Dep_Dates = [2454129.5:2454239.5];
Arr_Dates = [2456917.5:2457517.5];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Jupiter,V_Jupiter] = calcposvel(EphemD.Jupiter.a,EphemD.Jupiter.e,EphemD.Jupiter.i,EphemD.Jupiter.Omega,EphemD.Jupiter.w,EphemD.Jupiter.nu,mu);
        [R_Pluto,V_Pluto] = calcposvel(EphemA.Pluto.a,EphemA.Pluto.e,EphemA.Pluto.i,EphemA.Pluto.Omega,EphemA.Pluto.w,EphemA.Pluto.nu,mu);
        [V_initial,V_final] = solvelambert(R_Jupiter,R_Pluto,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Jupiter(i,j) = norm(V_initial - V_Jupiter);
        vinf_Pluto(i,j) = norm(V_final - V_Pluto);
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

V_inf_J_contours = [16,16.3,16.7,16.9,17,17.2,17.5,17.8,18,18.2,18.5,18.7,18.9,19,19.3,19.5,19.8,20,20.3,20.5,20.7,21];
V_inf_P_contours = [10,11,11.3,11.6,11.8,12,12.3,12.6,12.8,13,13.3,13.5,13.8,14.2,14.5,15,15.2,15.6,16];

dt_contours = [260,280,300,310,335,350,370,400,420,450,480,510,520,530,550,580];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,vinf_Jupiter,V_inf_J_contours,'r');
clabel(cs1,h1)
[cs2,h2] = contour(x_vals,y_vals,vinf_Pluto,V_inf_P_contours,'b');
clabel(cs2,h2)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k');
clabel(cs3,h3)
plot(30,300,'g*','MarkerSize',15)
xlabel('Days Past 29 Jan 2007')
ylabel('Days Past 17 Sept 2014')
legend('$V_{\infty,Jupiter}$,$\frac{km}{s}$','$V_{\infty,Pluto}$,$\frac{km}{s}$','TOF, days','NH Traj. Dates','Interpreter','latex')
title('Jupiter-Pluto New Horizons Pork Chop Plot')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Problem 6 -- Build search function that identifies candidate trajectories
% with valid gravity swingbys
clear C3 vinf_Jupiter dt 
mu = 1.32712440018e11; %mu of sun
Dep_Date = 2453745.29167;
Arr_Dates = [2454129.5:.5:2454239.5];
P_Arr_Dates = [2456917.5:.5:2457517.5];
JupRad = 71492;
trajcnt = 1;
validC3 = [];
trackj = [];

for i = 1:numel(Arr_Dates)
    EphemD = Meeus(Dep_Date);
    EphemA = Meeus(Arr_Dates(i));
    [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
    [R_Jupiter,V_Jupiter] = calcposvel(EphemA.Jupiter.a,EphemA.Jupiter.e,EphemA.Jupiter.i,EphemA.Jupiter.Omega,EphemA.Jupiter.w,EphemA.Jupiter.nu,mu);
    [V_initial,V_final] = solvelambert(R_Earth,R_Jupiter,(Arr_Dates(i) - Dep_Date)*86400,0);
    vinf_Jupiter(i) = norm(V_final - V_Jupiter);
    C3(i) = (norm(V_initial - V_Earth))^2;
    dt(i) = Arr_Dates(i) - Dep_Date;
    if C3(i) <= 180
        Dep_Jupiter = Arr_Dates(i); %Use viable Jupiter arrival date as Jupiter departure date
        for j = 1:numel(P_Arr_Dates)
            EphemD_J = Meeus(Dep_Jupiter);
            EphemA_P = Meeus(P_Arr_Dates(j));
            [R_Jupiter,V_Jupiter] = calcposvel(EphemD_J.Jupiter.a,EphemD_J.Jupiter.e,EphemD_J.Jupiter.i,EphemD_J.Jupiter.Omega,EphemD_J.Jupiter.w,EphemD_J.Jupiter.nu,mu);
            [R_Pluto,V_Pluto] = calcposvel(EphemA_P.Pluto.a,EphemA_P.Pluto.e,EphemA_P.Pluto.i,EphemA_P.Pluto.Omega,EphemA_P.Pluto.w,EphemA_P.Pluto.nu,mu);
            [V_initial,V_final] = solvelambert(R_Jupiter,R_Pluto,(P_Arr_Dates(j) - Dep_Jupiter)*86400,0);
            vinf_Jupiter_out = norm(V_initial - V_Jupiter);
            vinf_Pluto_in = norm(V_final - V_Pluto);
            %check if incoming and outgoing excess velocities are close
            %enough and if pluto excess velocity is low enough
            if abs(norm(vinf_Jupiter_out) - norm(vinf_Jupiter(i))) < 0.05
                if norm(vinf_Pluto_in) <= 14.5
                    num = dot(vinf_Jupiter(i),vinf_Jupiter_out);
                    denom = norm(vinf_Jupiter(i))*norm(vinf_Jupiter_out);
                    Psi = acos(num/denom);

                    rp = (mu_J/norm(vinf_Jupiter(i))^2)*((1/(cos((pi-Psi)/2))) - 1);

                    flyby_alt = rp - JupRad;

                    if flyby_alt >= 2144760
                        %trajectory meets all criteria, save params.
                        VTP.Jupiter_Arrival(trajcnt) = Arr_Dates(i);
                        VTP.Pluto_Arrival(trajcnt) = P_Arr_Dates(j);
                        VTP.C3(trajcnt) = C3(i);
                        VTP.V_inf_J_in(trajcnt,:) = vinf_Jupiter(i);
                        VTP.V_inf_J_out(trajcnt,:) = vinf_Jupiter_out;
                        VTP.V_inf_P_in(trajcnt,:) = vinf_Pluto_in;
                        VTP.flyby_alt(trajcnt) = flyby_alt;
                        VTP.flyby_Vinf_diff(trajcnt) = abs(VTP.V_inf_J_out(trajcnt) - VTP.V_inf_J_in(trajcnt));
                        trackj(trajcnt) = j;
                        trajcnt = trajcnt + 1;
                    end
                end
            end
        end
    end
end    
