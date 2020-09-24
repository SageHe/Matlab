clear all;close all;clc
Data1 = load('AM_Unit4_MOI_6mNm_G7');
Data1([1:107],:) = [];
time = Data1(:,1)./1000; %seconds
com_T = Data1(:,2); %commanded torque
w = abs(Data1(:,3))*((2*pi)/(60)); %angular velocity
current = Data1(:,4); %current through motor
%current conversion - 25.5 mNm/A
A_conv = 25.5;
%multiply current by oconversion to get actual torque on motor
app_T = (current*A_conv)./1000; %torque in Nm
alpha = gradient(time,w);
I1 = app_T./alpha;

% Second Data Set
Data2 = load('AM_Unit4_MOI_8mNm_G7');
Data2([1:109],:) = [];
time = Data2(:,1)./1000; %seconds
com_T = Data2(:,2); %commanded torque
w = abs(Data2(:,3))*((2*pi)/(60)); %angular velocity
current = Data2(:,4); %current through motor
%current conversion - 25.5 mNm/A
A_conv = 25.5;
%multiply current by oconversion to get actual torque on motor
app_T = (current*A_conv)./1000; %torque in mNm
alpha = gradient(time,w);
I2 = app_T./alpha;

% Third data set
Data3 = load('AM_Unit4_MOI_9mNm_G7');
Data3([1:104],:) = [];
time = Data3(:,1)./1000; %seconds
com_T = Data3(:,2); %commanded torque
w = abs(Data3(:,3))*((2*pi)/(60)); %angular velocity
current = Data3(:,4); %current through motor
%current conversion - 25.5 mNm/A
A_conv = 25.5;
%multiply current by oconversion to get actual torque on motor
app_T = (current*A_conv)./1000; %torque in mNm
alpha = gradient(time,w);
I3 = app_T./alpha;

I = mean([mean(I1) mean(I2) mean(I3)]);
