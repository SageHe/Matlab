clear all;close all;clc
Data1 = xlsread('VelocityVoltage_S011_G01.xlsx');

Atmos_pres1 = Data1(:,1); %Atmospheric Pressure in Pascals
Atmos_temp1 = Data1(:,2); %Atmospheric temperature in Kelvin
Airspeed_diff_pres1 = Data1(:,3); %Airspeed differential pressure in Pascals
Aux_diff_pres1 = Data1(:,4); %Aux differential pressure in Pascals
ELD_probe_x1 = Data1(:,5); %ELD probe x axis in mm
ELD_probe_y1 = Data1(:,6); %ELD probe y axis in mm
Voltage1 = Data1(:,7); %Voltge from measurements

Area_ratio = 9.5; %Given in lab document

V1_infinity_1 = Pitot_Static_EQ(Airspeed_diff_pres1,Atmos_temp1,Atmos_pres1);
% V1_infinity_2 = Venturi_Tube_Eq(Aux_diff_pres1,Atmos_temp1,Atmos_pres1,Area_ratio);

Data2 = xlsread('VelocityVoltage_S012_G03.xlsx');

Atmos_pres2 = Data2(:,1); %Atmospheric Pressure in Pascals
Atmos_temp2 = Data2(:,2); %Atmospheric temperature in Kelvin
Airspeed_diff_pres2 = Data2(:,3); %Airspeed differential pressure in Pascals
Aux_diff_pres2 = Data2(:,4); %Aux differential pressure in Pascals
ELD_probe_x2 = Data2(:,5); %ELD probe x axis in mm
ELD_probe_y2 = Data2(:,6); %ELD probe y axis in mm
Voltage2 = Data2(:,7); %Voltge from measurements

Area_ratio = 9.5; %Given in lab document

V2_infinity_1 = Pitot_Static_EQ(Airspeed_diff_pres2,Atmos_temp2,Atmos_pres2);
% V2_infinity_2 = Venturi_Tube_Eq(Aux_diff_pres2,Atmos_temp2,Atmos_pres2);


% plot(Voltage1,V1_infinity_1)
% Vvec = linspace(min(Voltage1),max(Voltage1),101);
% coeff = polyfit(Voltage1, V1_infinity_1,1);
% V_inf = coeff(1)*Vvec + coeff(2);
% hold on
% plot(Vvec,V_inf)
% hold off


