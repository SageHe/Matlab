clear all;close all;clc
%% Part 1
% Construct PCPs for each segment of the mission to narrow down valid
% critical event dates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Launch-Venus PCP

%Define specified event dates
Event.LaunchJD = 2460744.5;
% Event.VGAJD = 2460866.5;
% Event.VGAJD = 2460870.5;
Event.VGAJD = 2460901.5;
% Event.EGA1JD = 2461140.5;
% Event.EGA1JD = 2461221.5;
Event.EGA1JD = 2461252.5;
% Event.EGA2JD = 2461840.5;
Event.EGA2JD = 2462348.226567;
% Event.SOIJD = 2463901.5;
Event.SOIJD = 2463953.5;

mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.LaunchJD - 68:Event.LaunchJD + 45];
Arr_Dates = [Event.VGAJD - 45:Event.VGAJD + 45];
LaunchDates = Dep_Dates;
V_Arr_Dates = Arr_Dates;
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
fprintf('done\n')
% [I,J] = find(vinf_Venus == 5.91399779657392);
% [I,J] = find(C3 == min(min(C3)));

% clear vinf_Venus C3
% Dep_Date = Dep_Dates(49);
% Arr_Date = Arr_Dates(30);
%         EphemD = Meeus(Dep_Date);
%         EphemA = Meeus(Arr_Date);
%         [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
%         [R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
%         [V_initial,V_final] = solvelambert(R_Earth,R_Venus,(Arr_Date - Dep_Date)*86400,0);
%         vinf_Venus = norm(V_final - V_Venus);
%         C3 = (norm(V_initial - V_Earth))^2;
%         dt = Arr_Date - Dep_Date;
% 

[I,J] = find(C3 == min(min(C3)));

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [6 7 8 9 13 15 18 22 25];
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
plot(55,55,'*','MarkerSize',30)
set(gca,'FontSize',20)
xlabel('Days Past 1 Jan. 2025')
ylabel('Days Past 30 June 2025')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Venus}$,$\frac{km}{s}$','TOF, days','Chosen Event Date','FontSize',18,'Interpreter','latex')
title('Earth-Venus Pork Chop Plot','FontSize',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Venus-EGA1 PCP
clear vinf_Venus

mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.VGAJD - 140:Event.VGAJD + 150];
Arr_Dates = [Event.EGA1JD - 140:Event.EGA1JD + 150];
V_Dep_Dates = Dep_Dates;
EGA1_Arr_Dates = Arr_Dates;
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
fprintf('done\n')
% clear vinf_Venus clear vinf_Earth
% Dep_Date = Dep_Dates(144);
% Arr_Date = Arr_Dates(125);
%         EphemD = Meeus(Dep_Date);
%         EphemA = Meeus(Arr_Date);
%         [R_Venus,V_Venus] = calcposvel(EphemD.Venus.a,EphemD.Venus.e,EphemD.Venus.i,EphemD.Venus.Omega,EphemD.Venus.w,EphemD.Venus.nu,mu);
%         [R_Earth,V_Earth] = calcposvel(EphemA.Earth.a,EphemA.Earth.e,EphemA.Earth.i,EphemA.Earth.Omega,EphemA.Earth.w,EphemA.Earth.nu,mu);
%         [V_initial,V_final] = solvelambert(R_Venus,R_Earth,(Arr_Date - Dep_Date)*86400,0);
%         vinf_Venus = norm(V_initial - V_Venus);
%         vinf_Earth = norm(V_final - V_Earth);
%         dt = Arr_Date - Dep_Date;


x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

V_inf_Venus_contours = [2.5 2.8 3 3.5 4 4.5 5 6 7 8 9 11 14];
V_inf_Earth_contours = [2.5 2.8 3 3.5 4 4.5 5 6 6.5 7 8 10 11 12 13 15];
dt_contours = [200 220 240 260 280 300 320 340 360 380 400 420 440 460 480 500 520 540 560];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,vinf_Venus,V_inf_Venus_contours,'r','LineWidth',1);
clabel(cs1,h1,'FontSize',16)
[cs2,h2] = contour(x_vals,y_vals,vinf_Earth,V_inf_Earth_contours,'b','LineWidth',1);
clabel(cs2,h2,'FontSize',16)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k','LineWidth',1);
clabel(cs3,h3,'FontSize',16)
plot(150,115,'*','MarkerSize',30)
set(gca,'FontSize',20)
xlabel('Days Past 27 Mar 2025')
ylabel('Days Past 13 Mar 2026')
legend('$V_{\infty,Venus}$,$\frac{km}{s}$','$V_{\infty,Earth}$,$\frac{km}{s}$','TOF, days','Chosen Event Date','FontSize',18,'Interpreter','latex')
title('Venus-EGA1 Pork Chop Plot','FontSize',24)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EGA2-SOI PCP
clear vinf_Earth dt
% Event.EGA2JD = 2461951.984378;
mu = 1.32712440018e11; %mu of sun
Dep_Dates = [Event.EGA2JD - 120:Event.EGA2JD + 120];
Arr_Dates = [Event.SOIJD - 300:Event.SOIJD + 300];
EGA2_Dep_Dates = Dep_Dates;
Sat_Arr_Dates = Arr_Dates;
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Saturn,V_Saturn] = calcposvel(EphemA.Saturn.a,EphemA.Saturn.e,EphemA.Saturn.i,EphemA.Saturn.Omega,EphemA.Saturn.w,EphemA.Saturn.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Saturn,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_Earth(i,j) = norm(V_initial - V_Earth);
        vinf_Saturn(i,j) = norm(V_final - V_Saturn);
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

V_inf_Earth_contours = [8.5 8.8 9 9.5 10 11 13 15];
V_inf_Saturn_contours = [5 5.3 5.5 5.8 6 6.3 6.8 7 7.5 8 9 10 12 15];
dt_contours = [1650 1680 1750 1800 1850 1900 1950 2000 2050 2100 2150 2200 2350 2400 2450 2500 2550 2600];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,vinf_Earth,V_inf_Earth_contours,'r','LineWidth',1);
clabel(cs1,h1,'FontSize',16)
[cs2,h2] = contour(x_vals,y_vals,vinf_Saturn,V_inf_Saturn_contours,'b','LineWidth',1);
clabel(cs2,h2,'FontSize',16)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k','LineWidth',1);
clabel(cs3,h3,'FontSize',16)
plot(95,465,'*','MarkerSize',30)
set(gca,'FontSize',20)
xlabel('Days Past 1 Apr 2029')
ylabel('Days Past 25 Feb 2033')
legend('$V_{\infty,Earth}$,$\frac{km}{s}$','$V_{\infty,Saturn}$,$\frac{km}{s}$','TOF, days','Chosen Event Date','FontSize',18,'Interpreter','latex')
title('EGA2-SOI Pork Chop Plot','FontSize',24)
%% use searching function to determine valid dates of flybys
clear C3 dt 
mu = 1.32712440018e11; %mu of sun
mu_S = mu;
LaunchDate = 2460688.5;
% Arr_Dates = [2454129.5:.5:2454239.5];
% P_Arr_Dates = [2456917.5:.5:2457517.5];
% JupRad = 71492;
SatRad = 60268;
trajcnt = 1;
validC3 = [];
trackj = [];

for i = 1:numel(V_Arr_Dates)
    for n = 1:numel(LaunchDates)
        EphemD = Meeus(LaunchDates(n));
        EphemA = Meeus(V_Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Venus,(V_Arr_Dates(i) - LaunchDates(n))*86400,0);
        vinf_Venus_in(i) = norm(V_final - V_Venus);
        C3(i) = (norm(V_initial - V_Earth))^2;
        dt(i) = V_Arr_Dates(i) - LaunchDates(n);
        if C3(i) <= 25
            Dep_Venus = V_Arr_Dates(i); %Use viable Venus arrival date as Venus departure date
            for j = 1:numel(EGA1_Arr_Dates)
                EphemD_V = Meeus(Dep_Venus);
                EphemA_E = Meeus(EGA1_Arr_Dates(j));
                [R_Venus,V_Venus] = calcposvel(EphemD_V.Venus.a,EphemD_V.Venus.e,EphemD_V.Venus.i,EphemD_V.Venus.Omega,EphemD_V.Venus.w,EphemD_V.Venus.nu,mu);
                [R_Earth,V_Earth] = calcposvel(EphemA_E.Earth.a,EphemA_E.Earth.e,EphemA_E.Earth.i,EphemA_E.Earth.Omega,EphemA_E.Earth.w,EphemA_E.Earth.nu,mu);
                [V_initial,V_final] = solvelambert(R_Venus,R_Earth,(EGA1_Arr_Dates(j) - Dep_Venus)*86400,0);
                vinf_Venus_out = norm(V_initial - V_Venus);
                vinf_EGA1_in = V_final - V_Earth;
                R_EGA1 = R_Earth;
                V_Earth_EGA1 = V_Earth;
                %check if incoming and outgoing excess velocities are close
                %enough 
                if abs(norm(vinf_Venus_out) - norm(vinf_Venus_in(i))) < 0.05
                    %if dv_inf is low enough, add to current EGA1 arrival to
                    %get to EGA2 arrival with resonant orbit
                    EGA2_Arr_Date = EGA1_Arr_Dates(j) + 3*365.242189;
                    for k = 1:numel(Sat_Arr_Dates)
                        EphemD = Meeus(EGA2_Arr_Date);
                        EphemA = Meeus(Sat_Arr_Dates(k));
                        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
                        [R_Saturn,V_Saturn] = calcposvel(EphemA.Saturn.a,EphemA.Saturn.e,EphemA.Saturn.i,EphemA.Saturn.Omega,EphemA.Saturn.w,EphemA.Saturn.nu,mu);
                        [V_initial,V_final] = solvelambert(R_Earth,R_Saturn,(Sat_Arr_Dates(k) - EGA2_Arr_Date)*86400,0);
                        vinf_EGA2_out = V_initial - V_Earth;
                        vinf_Saturn_in = norm(V_final - V_Saturn);
                        dt = Sat_Arr_Dates(i) - EGA2_Arr_Date;
                        R_EGA2 = R_Earth;
                        V_Earth_EGA2 = V_Earth;

                    if norm(vinf_Saturn_in) < 9
                        if abs(norm(vinf_EGA2_out) - norm(vinf_EGA1_in)) < 0.03
                        %dont't need these requirements, consider replacing
                        %with res. orbit requirements over overlaying rp
                        %values?
                            mu_E = 3.98600433e5;
                            mu_V = 3.24858599e5;
                            AU = 1.49597870700e8;
                            mu_Sat = 3.7931208e7;

                            Re = 6378.14;
                            Rv = 6051.8;
                            Rs = 60268;

                            %Construct periapse radius of each gravity assist as func. of phi
                            %first flyby
                            P_Earth = 3*365.242189*86400;
                            a = ((P_Earth/(2*pi))^2*mu_S)^(1/3);
                            V_sc_post_EGA1 = sqrt(mu_S*((2/norm(R_EGA1)) - (1/a)));
                            V_inf = (norm(vinf_EGA2_out) + norm(vinf_EGA1_in))/2;
                            V_sc_sun = sqrt((2*mu_S)/norm(R_EGA1) - (mu_S/a));

                            num = -norm(V_sc_sun)^2 + V_inf^2 + norm(V_Earth_EGA1)^2;
                            denom = 2*V_inf*norm(V_Earth_EGA1);
                            theta = acos((num/denom));
                            %construct VNC frame transformation
                            V_hat = V_Earth_EGA1/norm(V_Earth_EGA1);
                            N_hat = (cross(R_EGA1,V_Earth_EGA1))/norm(cross(R_EGA1,V_Earth_EGA1));
                            C_hat = cross(V_hat,N_hat);
                            T = [V_hat N_hat C_hat];

                            phi = [0:0.01:2*pi];

                            for m = 1:numel(phi)
                                V_inf_out_EGA1 = V_inf*[cos(pi - theta);sin(pi - theta)*cos(phi(m));-sin(pi - theta)*sin(phi(m))];
                                V_inf_out_EGA1(:,m) = T*V_inf_out_EGA1;
                                V_inf_in_EGA2(:,m) = V_inf_out_EGA1(:,m) + V_Earth_EGA1 - V_Earth_EGA2;
                                Psi_EGA1(m) = acos((dot(vinf_EGA1_in,V_inf_out_EGA1(:,m)))/(norm(vinf_EGA1_in)*norm(V_inf_out_EGA1(:,m))));
                                Psi_EGA2(m) = acos((dot(V_inf_in_EGA2(:,m),vinf_EGA2_out))/(norm(V_inf_in_EGA2(:,m))*norm(vinf_EGA2_out)));
                                rp_EGA1(m) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA1(m))/2) - 1));
                                rp_EGA2(m) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA2(m))/2) - 1));
                            end
                            phiEGA1 = find(rp_EGA1 > (Re+300));
                            phiEGA2 = find(rp_EGA2 > (Re+300));
                            acceptable_phi = intersect(phiEGA1,phiEGA2);
                            plot_acc_phix = [rad2deg(phi(min(acceptable_phi))) rad2deg(phi(min(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) ];
                            plot_acc_phiy = [0 12000 12000 0];    
                            rplim = (Re+300)*ones(numel(phi),1);
                            
                            if ~isempty(acceptable_phi)
%                                 fprintf('good traj.\n')
%                                 figure
%                                 hold on
%                                 plot(phi*(180/pi),rp_EGA1,'LineWidth',1)
%                                 plot(phi*(180/pi),rp_EGA2,'LineWidth',1)
%                                 plot(phi*(180/pi),rplim,'--k','LineWidth',1)
%                                 patch(plot_acc_phix,plot_acc_phiy,'gree','FaceAlpha',0.2)
%                                 set(gca,'FontSize',20)
%                                 xlim([0 360])
%                                 ylim([0 12000])
%                                 grid on
%                                 grid minor
% 
%                                 xlabel('\phi (Degrees)')
%                                 ylabel('r_p (km)')
%                                 title('Periapse Radius r_p VS \phi','FontSize',24)
%                                 legend('EGA1 r_p','EGA2 r_p','Minimum r_p','Acceptable \phi region','FontSize',16)


%                             figure
%                             hold on
%                             plot(phi*(180/pi),rp_EGA1,'LineWidth',1)
%                             plot(phi*(180/pi),rp_EGA2,'LineWidth',1)
%                             plot(phi*(180/pi),rplim,'--k','LineWidth',1)
%     %                         patch(plot_acc_phix,plot_acc_phiy,'gree','FaceAlpha',0.2)
%     %                         set(gca,'FontSize',20)
%     %                         xlim([0 360])
%     %                         ylim([0 12000])
%     %                         grid on
%     %                         grid minor
%     % 
%     %                         xlabel('\phi (Degrees)')
%     %                         ylabel('r_p (km)')
%     %                         title('Periapse Radius r_p VS \phi','FontSize',24)
%     %                         legend('EGA1 r_p','EGA2 r_p','Minimum r_p','Acceptable \phi region','FontSize',16)
    % 
                                VTP.LaunchDateJD(trajcnt) = LaunchDates(n);
                                VTP.VGA(trajcnt) = V_Arr_Dates(i);
                                VTP.EGA1Date(trajcnt) = EGA1_Arr_Dates(j);
                                VTP.EGA2Date(trajcnt) = EGA2_Arr_Date;
                                VTP.SaturnDate(trajcnt) = Sat_Arr_Dates(k);
                                VTP.C3(trajcnt) = C3(i);
                                VTP.TOF(trajcnt) = Sat_Arr_Dates(k) - LaunchDates(n);
                                VTP.vinf_V_in(trajcnt) = vinf_Venus_in(i);
                                VTP.vinf_V_out(trajcnt) = vinf_Venus_out;
                                VTP.vinf_EGA1_in(trajcnt) = norm(vinf_EGA1_in);
                                VTP.vinf_EGA2_out(trajcnt) = norm(vinf_EGA2_out);
                                VTP.vinf_Saturn_in(trajcnt) = vinf_Saturn_in;
                                VTP.res_phi{trajcnt} = acceptable_phi;
                                VTP.res_rp_EGA1(trajcnt,:) = rp_EGA1;
                                VTP.res_rp_EGA2(trajcnt,:) = rp_EGA2;
                                trajcnt = trajcnt + 1;
                            end
                        end
                    end
                    end
                end
            end
        end
    end
        i
end    
fprintf('Done searching, %d valid trajectories found\n',trajcnt)
%%
V_Arr_Dates = VTP.VGA(63);
LaunchDates = VTP.LaunchDateJD(63);
EGA1_Arr_Dates = VTP.EGA1Date(63);
EGA2_Arr_Dates = VTP.EGA2Date(63);
Sat_Arr_Dates = VTP.SaturnDate(63);

for i = 1:numel(V_Arr_Dates)
    for n = 1:numel(LaunchDates)
        EphemD = Meeus(LaunchDates(n));
        EphemA = Meeus(V_Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Venus,V_Venus] = calcposvel(EphemA.Venus.a,EphemA.Venus.e,EphemA.Venus.i,EphemA.Venus.Omega,EphemA.Venus.w,EphemA.Venus.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Venus,(V_Arr_Dates(i) - LaunchDates(n))*86400,0);
        vinf_Venus_in(i) = norm(V_final - V_Venus);
        C3(i) = (norm(V_initial - V_Earth))^2;
        dt(i) = V_Arr_Dates(i) - LaunchDates(n);
        if C3(i) <= 25
            Dep_Venus = V_Arr_Dates(i); %Use viable Venus arrival date as Venus departure date
            for j = 1:numel(EGA1_Arr_Dates)
                EphemD_V = Meeus(Dep_Venus);
                EphemA_E = Meeus(EGA1_Arr_Dates(j));
                [R_Venus,V_Venus] = calcposvel(EphemD_V.Venus.a,EphemD_V.Venus.e,EphemD_V.Venus.i,EphemD_V.Venus.Omega,EphemD_V.Venus.w,EphemD_V.Venus.nu,mu);
                [R_Earth,V_Earth] = calcposvel(EphemA_E.Earth.a,EphemA_E.Earth.e,EphemA_E.Earth.i,EphemA_E.Earth.Omega,EphemA_E.Earth.w,EphemA_E.Earth.nu,mu);
                [V_initial,V_final] = solvelambert(R_Venus,R_Earth,(EGA1_Arr_Dates(j) - Dep_Venus)*86400,0);
                vinf_Venus_out = norm(V_initial - V_Venus);
                vinf_EGA1_in = V_final - V_Earth;
                R_EGA1 = R_Earth;
                V_Earth_EGA1 = V_Earth;
                %check if incoming and outgoing excess velocities are close
                %enough 
                if abs(norm(vinf_Venus_out) - norm(vinf_Venus_in(i))) < 0.05
                    %if dv_inf is low enough, add to current EGA1 arrival to
                    %get to EGA2 arrival with resonant orbit
                    EGA2_Arr_Date = EGA1_Arr_Dates(j) + 3*365.242189;
                    for k = 1:numel(Sat_Arr_Dates)
                        EphemD = Meeus(EGA2_Arr_Date);
                        EphemA = Meeus(Sat_Arr_Dates(k));
                        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
                        [R_Saturn,V_Saturn] = calcposvel(EphemA.Saturn.a,EphemA.Saturn.e,EphemA.Saturn.i,EphemA.Saturn.Omega,EphemA.Saturn.w,EphemA.Saturn.nu,mu);
                        [V_initial,V_final] = solvelambert(R_Earth,R_Saturn,(Sat_Arr_Dates(k) - EGA2_Arr_Date)*86400,0);
                        vinf_EGA2_out = V_initial - V_Earth;
                        vinf_Saturn_in = norm(V_final - V_Saturn);
                        dt = Sat_Arr_Dates(i) - EGA2_Arr_Date;
                        R_EGA2 = R_Earth;
                        V_Earth_EGA2 = V_Earth;

                    if norm(vinf_Saturn_in) < 9
                        if abs(norm(vinf_EGA2_out) - norm(vinf_EGA1_in)) < 0.03
                        %dont't need these requirements, consider replacing
                        %with res. orbit requirements over overlaying rp
                        %values?
                            mu_E = 3.98600433e5;
                            mu_V = 3.24858599e5;
                            AU = 1.49597870700e8;
                            mu_Sat = 3.7931208e7;

                            Re = 6378.14;
                            Rv = 6051.8;
                            Rs = 60268;

                            %Construct periapse radius of each gravity assist as func. of phi
                            %first flyby
                            P_Earth = 3*365.242189*86400;
                            a = ((P_Earth/(2*pi))^2*mu_S)^(1/3);
                            V_sc_post_EGA1 = sqrt(mu_S*((2/norm(R_EGA1)) - (1/a)));
                            V_inf = (norm(vinf_EGA2_out) + norm(vinf_EGA1_in))/2;
                            V_sc_sun = sqrt((2*mu_S)/norm(R_EGA1) - (mu_S/a));

                            num = -norm(V_sc_sun)^2 + V_inf^2 + norm(V_Earth_EGA1)^2;
                            denom = 2*V_inf*norm(V_Earth_EGA1);
                            theta = acos((num/denom));
                            %construct VNC frame transformation
                            V_hat = V_Earth_EGA1/norm(V_Earth_EGA1);
                            N_hat = (cross(R_EGA1,V_Earth_EGA1))/norm(cross(R_EGA1,V_Earth_EGA1));
                            C_hat = cross(V_hat,N_hat);
                            T = [V_hat N_hat C_hat];

                            phi = [0:0.01:2*pi];

                            for m = 1:numel(phi)
                                V_inf_out_EGA1 = V_inf*[cos(pi - theta);sin(pi - theta)*cos(phi(m));-sin(pi - theta)*sin(phi(m))];
                                V_inf_out_EGA1(:,m) = T*V_inf_out_EGA1;
                                V_inf_in_EGA2(:,m) = V_inf_out_EGA1(:,m) + V_Earth_EGA1 - V_Earth_EGA2;
                                Psi_EGA1(m) = acos((dot(vinf_EGA1_in,V_inf_out_EGA1(:,m)))/(norm(vinf_EGA1_in)*norm(V_inf_out_EGA1(:,m))));
                                Psi_EGA2(m) = acos((dot(V_inf_in_EGA2(:,m),vinf_EGA2_out))/(norm(V_inf_in_EGA2(:,m))*norm(vinf_EGA2_out)));
                                rp_EGA1(m) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA1(m))/2) - 1));
                                rp_EGA2(m) = (mu_E/V_inf^2)*((1/cos((pi - Psi_EGA2(m))/2) - 1));
                            end
                            phiEGA1 = find(rp_EGA1 > (Re+300));
                            phiEGA2 = find(rp_EGA2 > (Re+300));
                            acceptable_phi = intersect(phiEGA1,phiEGA2);
                            plot_acc_phix = [rad2deg(phi(min(acceptable_phi))) rad2deg(phi(min(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) rad2deg(phi(max(acceptable_phi))) ];
                            plot_acc_phiy = [0 12000 12000 0];    
                            rplim = (Re+300)*ones(numel(phi),1);
                            
                            if ~isempty(acceptable_phi)
                                fprintf('good traj.\n')
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


%                             figure
%                             hold on
%                             plot(phi*(180/pi),rp_EGA1,'LineWidth',1)
%                             plot(phi*(180/pi),rp_EGA2,'LineWidth',1)
%                             plot(phi*(180/pi),rplim,'--k','LineWidth',1)
%     %                         patch(plot_acc_phix,plot_acc_phiy,'gree','FaceAlpha',0.2)
%     %                         set(gca,'FontSize',20)
%     %                         xlim([0 360])
%     %                         ylim([0 12000])
%     %                         grid on
%     %                         grid minor
%     % 
%     %                         xlabel('\phi (Degrees)')
%     %                         ylabel('r_p (km)')
%     %                         title('Periapse Radius r_p VS \phi','FontSize',24)
%     %                         legend('EGA1 r_p','EGA2 r_p','Minimum r_p','Acceptable \phi region','FontSize',16)
    % 
                                VTP.LaunchDateJD(trajcnt) = LaunchDates(n);
                                VTP.VGA(trajcnt) = V_Arr_Dates(i);
                                VTP.EGA1Date(trajcnt) = EGA1_Arr_Dates(j);
                                VTP.EGA2Date(trajcnt) = EGA2_Arr_Date;
                                VTP.SaturnDate(trajcnt) = Sat_Arr_Dates(k);
                                VTP.C3(trajcnt) = C3(i);
                                VTP.TOF(trajcnt) = Sat_Arr_Dates(k) - LaunchDates(n);
                                VTP.vinf_V_in(trajcnt) = vinf_Venus_in(i);
                                VTP.vinf_V_out(trajcnt) = vinf_Venus_out;
                                VTP.vinf_EGA1_in(trajcnt) = norm(vinf_EGA1_in);
                                VTP.vinf_EGA2_out(trajcnt) = norm(vinf_EGA2_out);
                                VTP.vinf_Saturn_in(trajcnt) = vinf_Saturn_in;
                                VTP.res_phi{trajcnt} = acceptable_phi;
                                VTP.res_rp_EGA1(trajcnt,:) = rp_EGA1;
                                VTP.res_rp_EGA2(trajcnt,:) = rp_EGA2;
                                trajcnt = trajcnt + 1;
                            end
                        end
                    end
                    end
                end
            end
        end
    end
        i
end    

% fprintf('Searching Function Done\n')
%% Part 2 
% Construct Earth-Earth resonant orbit plot and state phi value that was
% selected and why

mu_S = 1.32712440018e11;
mu_E = 3.98600433e5;
mu_V = 3.24858599e5;
AU = 1.49597870700e8;
mu_Sat = 3.7931208e7;

Re = 6378.14;
Rv = 6051.8;
Rs = 60268;

% Event.Launch = 2458951.5;
% Event.VGA1 = 2459121.5;
% Event.EGA1 = 2459438.5;
% Event.EGA2 = 2460534.226567;
% Event.JOI = 2461300.75;

% Event.VGAJD = 2460870.5;
% Event.EGA1JD = 2461221.5;
% Event.EGA2JD = 2462317.226567;

% Calculate hyperbolic excess velocity of S/C upon arrivalt at EGA1
VGA1Ephem = Meeus(Event.VGAJD);
EGA1Ephem = Meeus(Event.EGA1JD);

[R_VGA1,Vp_VGA1] = calcposvel(VGA1Ephem.Venus.a,VGA1Ephem.Venus.e,VGA1Ephem.Venus.i,VGA1Ephem.Venus.Omega,VGA1Ephem.Venus.w,VGA1Ephem.Venus.nu,mu_S);
[R_EGA1,Vp_EGA1] = calcposvel(EGA1Ephem.Earth.a,EGA1Ephem.Earth.e,EGA1Ephem.Earth.i,EGA1Ephem.Earth.Omega,EGA1Ephem.Earth.w,EGA1Ephem.Earth.nu,mu_S);
% V_inf_in = vf - v0V, V_inf_out = v0_V_FB - v0V;
[v_sc_VGA1,v_sc_EGA1] = solvelambert(R_VGA1,R_EGA1,(Event.EGA1JD - Event.VGAJD)*86400,0);

V_inf_in_EGA1 = v_sc_EGA1 - Vp_EGA1;

%Calculate hyperbolic excess velocity of S/C up departure at EGA2
EGA2Ephem = Meeus(Event.EGA2JD);
SOIEphem = Meeus(Event.SOIJD);
% 
[R_EGA2,Vp_EGA2] = calcposvel(EGA2Ephem.Earth.a,EGA2Ephem.Earth.e,EGA2Ephem.Earth.i,EGA2Ephem.Earth.Omega,EGA2Ephem.Earth.w,EGA2Ephem.Earth.nu,mu_S);
[R_SOI,Vp_SOI] = calcposvel(SOIEphem.Saturn.a,SOIEphem.Saturn.e,SOIEphem.Saturn.i,SOIEphem.Saturn.Omega,SOIEphem.Saturn.w,SOIEphem.Saturn.nu,mu_S);

[v_sc_EGA2,v_sc_SOI] = solvelambert(R_EGA2,R_SOI,(Event.SOIJD - Event.EGA2JD)*86400,0);

V_inf_out_EGA2 = v_sc_EGA2 - Vp_EGA2;

V_inf_in_SOI = v_sc_SOI - Vp_SOI;

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
