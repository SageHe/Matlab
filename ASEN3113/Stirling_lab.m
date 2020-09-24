%%Thermo lab 1-Stirling Engine
%Author - Sage Herrin
%Created - 9/13/18
%ID - 106071909

clear all;close all;clc
%Read in data from text file 
data = readtable('Group4_012_8.txt');
time = table2array(data(:,1));
pressure = table2array(data(:,2));
top_of_top = table2array(data(:,3));
bottom_of_top = table2array(data(:,4));
top_of_bottom = table2array(data(:,5));
bottom_of_bottom = table2array(data(:,6));
current = table2array(data(:,7));
opt_switch = table2array(data(:,8));
data = table2array(data);
%Separate different trials of engine
%Trial 1
temp1 = find(time < 6);
data1 = data(1:temp1(end),:);
time1 = (data1(:,1));
pressure1 = (data1(:,2));
top_of_top1 = (data1(:,3));
bottom_of_top1 = (data1(:,4));
top_of_bottom1 = (data1(:,5));
bottom_of_bottom1 = (data1(:,6));
current1 = (data1(:,7));
opt_switch1 = (data1(:,8));
%Trial 2
temp2 = find(time > 143);
data2 = data((9076:39600),:);
time2 = (data2(:,1));
pressure2 = (data2(:,2));
top_of_top2 = (data(:,3));
bottom_of_top2 = (data2(:,4));
top_of_bottom2 = (data2(:,5));
bottom_of_bottom2 = (data2(:,6));
current2 = (data2(:,7));
opt_switch2 = (data2(:,8));
%Trial 3
data3 = data((39601:end),:);
time3 = (data3(:,1));
pressure3 = (data3(:,2));
top_of_top3 = (data3(:,3));
bottom_of_top3 = (data3(:,4));
top_of_bottom3 = (data3(:,5));
bottom_of_bottom3 = (data3(:,6));
current3 = (data3(:,7));
opt_switch3 = (data3(:,8));

%RPM Calculations trial 1
dt = diff(opt_switch1);
dt = find(dt ~= 0);
RPM1 = [];
for i = 1:((length(dt)/2) - 1)
    RPM1 = [RPM1 time1(dt((2*i) + 1)) - time1(dt((2*i) - 1))];
end
RPM1 = 1./RPM1;
RPM1 = RPM1.*60;
%RPM Calculations trial 2
dt = diff(opt_switch2);
dt = find(dt ~= 0);
RPM2 = [];
for i = 1:((length(dt)/2) - 1)
    RPM2 = [RPM2 time2(dt((2*i) + 1)) - time2(dt((2*i) - 1))];
end
RPM2 = 1./RPM2;
RPM2 = RPM2.*60;
%RPM Calculations trial 3
dt = diff(opt_switch3);
dt = find(dt ~= 0);
RPM3 = [];
for i = 1:((length(dt)/2) - 1)
    RPM3 = [RPM3 time3(dt((2*i) + 1)) - time3(dt((2*i) - 1))];
end
RPM3 = 1./RPM3;
RPM3 = RPM3.*60;
%Average temps trial 1 
mean_bottom1 = mean([mean(bottom_of_bottom1) mean(top_of_bottom1)]);
mean_top1 = mean([mean(bottom_of_top1) mean(top_of_top1)]);
%Average temps trial 2
mean_bottom2 = mean([mean(bottom_of_bottom2) mean(top_of_bottom2)]);
mean_top2 = mean([mean(bottom_of_top2) mean(top_of_top2)]);
%Average temps trial 3
mean_bottom3 = mean([mean(bottom_of_bottom3) mean(top_of_bottom3)]);
mean_top3 = mean([mean(bottom_of_top3) mean(top_of_top3)]);
%Power piston diameter = 15mm
[NUM,TXT,RAW] = xlsread('Small_bottomface_lin_disp_trial1_87.46RPM.xlsx');
disp = NUM(:,3);
disp = abs(disp);
disp = (max(disp) - min(disp))/1000; %m
piston_diameter = .015; %m
piston_area = pi*(piston_diameter/2)^2;
working_vol_small = piston_area*disp;
%Big piston, large diameter = 144mm, small diameter = 140mm, foam thickness
%= 11mm, piston height = 21mm
working_vol_large = pi*(.072)^2*.021 - pi*(.07)^2*.011;
working_vol_tot = working_vol_large + working_vol_small;
%%Trial 1
PV1_time = time1(7364:8505);
PV1_time = PV1_time - PV1_time(1);
PV1_pres = pressure1(7364:8505);
PV1_pres = PV1_pres + 12.2;
PV1_pres = PV1_pres*6894.757; %Pascals

plot(PV1_time,PV1_pres)  

PP_data = xlsread('Small_bottomface_lin_disp_trial1_87.46RPM.xlsx');
figure
plot(PP_data(:,2),PP_data(:,3));
hold on;
y_line = linspace(PP_data(1,3),PP_data(1,3),2101);
plot(PP_data(:,2),y_line,'b');
relevant_data = [PP_data(185:277,2), PP_data(185:277,3)];
plot(relevant_data(:,1),relevant_data(:,2),'r');
 
time = relevant_data(:,1);
time = time-(time(1));

 
displace = relevant_data(:,2);
bottom = min(displace);
displace = displace - bottom;
displace = displace/1000;
displace = displace*(.015^2)*pi;
wv = working_vol_large + displace; 
figure
plot(time,wv);

t1 = linspace(1,length(PV1_pres),length(wv));
t2 = 1:length(PV1_pres);

PV_VOL1 = interp1(t1,wv,t2);

figure
plot(PV_VOL1,PV1_pres)
title('PV Diagram for 8 Degree Temperature Difference')
xlabel('Volume [m^3]')
ylabel('Pressure [Pa]')

%%Trial 2
PV2_time = time2(28698:29573);
PV2_time = PV2_time - PV2_time(1);
PV2_pres = pressure2(28698:29573);
PV2_pres = PV2_pres + 12.2;
PV2_pres = PV2_pres*6894.757; %Pascals

figure
plot(PV2_time,PV2_pres)  

figure
PP_data_2 = xlsread('small_bottomface_lin_disp_trial2_112.89RPM.xlsx');
plot(PP_data_2(:,2),PP_data_2(:,3));
hold on;
y_line = linspace(PP_data_2(1,3),PP_data_2(1,3),2101);
plot(PP_data_2(:,2),y_line,'b');
rd_2 = [PP_data_2(532:619,2), PP_data_2(532:619,3)];
plot(rd_2(:,1),rd_2(:,2));
 

t2 = rd_2(:,1);
t2 = t2-(t2(1));
 

disp2 = rd_2(:,2);
mini = min(disp2);
disp2 = disp2-mini;
disp2 = disp2/1000;
disp2 = disp2*(.015^2)*pi;
wv2 = working_vol_large + disp2;
figure
plot(t2,wv2)

t1_2 = linspace(1,length(PV2_pres),length(wv2));
t2_2 = 1:length(PV2_pres);

PV_VOL2 = interp1(t1_2,wv2,t2_2);
 
figure
plot(PV_VOL2,PV2_pres)
title('PV Diagram for 10 Degree Temperature Difference')
xlabel('Volume [m^3]')
ylabel('Pressure [Pa]')

%%Trial 3
PV3_time = time3(16547:17295);
PV3_time = PV3_time - PV3_time(1);
PV3_pres = pressure3(16547:17295);
PV3_pres = PV3_pres + 12.2;
PV3_pres = PV3_pres*6894.757; %Pascals

figure
plot(PV3_time,PV3_pres)  
figure
PP_data_3 = xlsread('small_bottomface_lin_disp_trial3_132.26RPM.xlsx');
plot(PP_data_3(:,2),PP_data_3(:,3));
hold on;
y_line = linspace(PP_data_3(1,3),PP_data_3(1,3),2101);
plot(PP_data_3(:,2),y_line,'b');
rd_3 = [PP_data_3(351:436,2), PP_data_3(351:436,3)];
plot(rd_3(:,1),rd_3(:,2));
 

t3 = rd_3(:,1);
t3 = t3-(t3(1));
 

disp3 = rd_3(:,2);
mini3 = min(disp3);
disp3 = disp3-mini3;
disp3 = disp3/1000;
disp3  = disp3*(.015^2 * pi);
wv3 = working_vol_large + disp3;
figure
plot(t3,wv3);

t1_3 = linspace(1,length(PV3_pres),length(wv3));
t2_3 = 1:length(PV3_pres);

PV_VOL3 = interp1(t1_3,wv3,t2_3);
 
figure
plot(PV_VOL3,PV3_pres)
title('PV Diagram for 12 Degree Temperature Difference')
xlabel('Volume [m^3]')
ylabel('Pressure [Pa]')

% plot(PV_VOL3(1),PV3_pres(1),'or')
%%Work integration 
%Trial 1
%left index
left = find(PV_VOL1 == min(PV_VOL1));
right = find(PV_VOL1 == max(PV_VOL1));
%lower left half
b_1 = trapz(PV_VOL1(1:left),PV1_pres(1:left))*-1;
%lower right half
c_1 = trapz(PV_VOL1(right:(end - 9)),PV1_pres(right:(end - 9)))*-1;
%upper half
a_1 = trapz(PV_VOL1(left:right),PV1_pres(left:right));
wnet_1 = a_1 - (c_1 + b_1);
%Trial 2
%left index
left = find(PV_VOL2 == min(PV_VOL2));
right = find(PV_VOL2 == max(PV_VOL2));
%lower left half
b_2 = trapz(PV_VOL2(1:left),PV2_pres(1:left))*-1;
%lower right half
c_2 = trapz(PV_VOL2(right:(end - 1)),PV2_pres(right:(end - 1)))*-1;
%upper half
a_2 = trapz(PV_VOL2(left:right),PV2_pres(left:right));
wnet_2 = a_2 - (c_2 + b_2);
%Trial 3
%left index
left = find(PV_VOL3 == min(PV_VOL3));
right = find(PV_VOL3 == max(PV_VOL3));
%lower left half
b_3 = trapz(PV_VOL3(1:left),PV3_pres(1:left))*-1;
%lower right half
c_3 = trapz(PV_VOL3(right:(end - 1)),PV3_pres(right:(end - 1)))*-1;
%upper half
a_3 = trapz(PV_VOL3(left:right),PV3_pres(left:right));
wnet_3 = a_3 - (c_3 + b_3);
%%Input work using P = IV
%Input work 1
I1 = mean(current1);
P1 = I1*5*(1/(mean(RPM1)/60)); %constant 5 volts
%Input work 2
I2 = mean(current2);
P2 = I2*5*(1/(mean(RPM2)/60));
%Input work 2
I3 = mean(current3);
P3 = I3*5*(1/(mean(RPM3)/60));
