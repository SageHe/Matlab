clear all;close all;clc
%% Lab 4 -- pork chop plot generation for Earth-Mars flight
%Question 2 -- generate pork chop plot for 2026 Earth/Mars opportunity.
%dep. dates range from JD=2461295 to 2461415 and arrival dates ranging from
%JD=2461530 to 2461710. Plot C3 at Earth, V_inf at Mars and TOF in days

mu = 1.32712440018e11; %mu of sun
Dep_Dates = [2461295:2461415];
Arr_Dates = [2461530:2461710];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Mars,V_Mars] = calcposvel(EphemA.Mars.a,EphemA.Mars.e,EphemA.Mars.i,EphemA.Mars.Omega,EphemA.Mars.w,EphemA.Mars.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Mars,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_mars(i,j) = norm(V_final - V_Mars);
        C3(i,j) = (norm(V_initial - V_Earth))^2;
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [9.5,10.5,12.5,15,20,35];
V_inf_contours = [-2,-1,0,1,1.5,2,2.5,3,3.5,4,5,7,8.5,11,12,12.5,13];
dt_contours = [200,225,235,265,280,295,310,330,350,415];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,C3,C3_contours,'r');
clabel(cs1,h1)
[cs2,h2] = contour(x_vals,y_vals,vinf_mars,V_inf_contours,'b');
clabel(cs2,h2)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k');
clabel(cs3,h3)
xlabel('Days Past 11 Sept 2026')
ylabel('Days Past 4 May 2027')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Mars}$,$\frac{km}{s}$','TOF, days','Interpreter','latex')
title('Earth-Mars 2026 Opportunity Pork Chop Plot')
%% Question 3 
%{
-- gen. pork chop plot for 2028 Earth to Mars Opp.. Use Earth dep. dates
from JD=2462075 to 246175 and Marss arrival dates from JD=2462280 to
2462460. Plot C3 @ Earth and V_inf at Mars and TOF in days.
%}
clearvars -except mu
Dep_Dates = [2462075:2462175];
Arr_Dates = [2462280:2462460];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Mars,V_Mars] = calcposvel(EphemA.Mars.a,EphemA.Mars.e,EphemA.Mars.i,EphemA.Mars.Omega,EphemA.Mars.w,EphemA.Mars.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Mars,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_mars(i,j) = norm(V_final - V_Mars);
        C3(i,j) = (norm(V_initial - V_Earth))^2;
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [9.5,10.5,12.5,15,20,35];
V_inf_contours = [-2,-1,0,1,1.5,2,2.5,3,3.5,4,5,7,8.5,11];
dt_contours = [200,220,235,265,280,295];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,C3,C3_contours,'r');
clabel(cs1,h1)
[cs2,h2] = contour(x_vals,y_vals,vinf_mars,V_inf_contours,'b');
clabel(cs2,h2)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k');
clabel(cs3,h3)
xlabel('Days Past 30 Oct 2028')
ylabel('Days Past 23 May 2029')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Mars}$,$\frac{km}{s}$','TOF, days','Interpreter','latex')
title('Earth-Mars 2028 Opportunity Pork Chop Plot')
%% Bonus Question
clearvars -except mu
Dep_Dates = [juliandate('01-Jul-2020','dd-mm-yyyy'):juliandate('30-Sep-2020','dd-mm-yyyy')];
Arr_Dates = [juliandate('30-jan-2021','dd-mm-yyyy'):juliandate('10-Mar-2021','dd-mm-yyyy')];
for j = 1:numel(Dep_Dates)
    for i = 1:numel(Arr_Dates)
        EphemD = Meeus(Dep_Dates(j));
        EphemA = Meeus(Arr_Dates(i));
        [R_Earth,V_Earth] = calcposvel(EphemD.Earth.a,EphemD.Earth.e,EphemD.Earth.i,EphemD.Earth.Omega,EphemD.Earth.w,EphemD.Earth.nu,mu);
        [R_Mars,V_Mars] = calcposvel(EphemA.Mars.a,EphemA.Mars.e,EphemA.Mars.i,EphemA.Mars.Omega,EphemA.Mars.w,EphemA.Mars.nu,mu);
        [V_initial,V_final] = solvelambert(R_Earth,R_Mars,(Arr_Dates(i) - Dep_Dates(j))*86400,0);
        vinf_mars(i,j) = norm(V_final - V_Mars);
        C3(i,j) = (norm(V_initial - V_Earth))^2;
        dt(i,j) = Arr_Dates(i) - Dep_Dates(j);
        i
        j
    end
end

x_vals = Dep_Dates - Dep_Dates(1);
y_vals = Arr_Dates - Arr_Dates(1);

C3_contours = [9.5,10.5,12.5,13,15,17,20,22,25];
V_inf_contours = [-2,-1,0,1,1.5,2,2.5,3,3.5,4,5,7,8.5,11];
dt_contours = [200,220,235,265,280,295];

figure
hold on
[cs1,h1] = contour(x_vals,y_vals,C3,C3_contours,'r');
clabel(cs1,h1)
[cs2,h2] = contour(x_vals,y_vals,vinf_mars,V_inf_contours,'b');
clabel(cs2,h2)
[cs3,h3] = contour(x_vals,y_vals,dt,dt_contours,'k');
clabel(cs3,h3)
xlabel('Days Past 1 Jan 2020')
ylabel('Days Past 28 Aug 2020')
legend('C3 $\frac{km^2}{s^2}$','$V_{\infty,Mars}$,$\frac{km}{s}$','TOF, days','Interpreter','latex')
title('Earth-Mars 2020 Opportunity Pork Chop Plot')