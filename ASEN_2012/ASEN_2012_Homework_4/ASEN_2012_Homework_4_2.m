clear all;close all;clc
%Sage Herrin, ASEN 2012, Due 10/26/17, 2PM lecture
%The purpose of this problem was to determine the water height in the tank
%by using the known area of the tank base along with the matlab functin
%'trapz' 

%Time in minutes from 0 to 10 and flow rate in cubic feet per minute given
%in table. Use matlab function trapz to determine the approximate water
%height at t = 10 minutes. 

%Use fundamental equations, i.e. volume of a cirular cylinder. Integrate
%flow rate over time interval to find total volume of water in tank and
%manipulate equation for volume of cylinder using known volume of water and
%given base area to solve for water height after time t = 10 minutes

%From table given in assignment

time = 0:10;
flow_rate = [0 80 130 150 150 160 165 170 160 140 120];

%The area of the tank base is given to bee 100 square feet, which can be
%manipulated to find the radiius of the tank base

Area = 100;
rad = sqrt(100 / pi);

%Once the radius is calculated, the equation for the volume of a circular
%cylinder can be manipulated to find the height of the water in the tank
%given the total volume of the tank, calculated using the built-in matlab
%function trapz

vol = trapz(time,flow_rate);

height = vol / (pi * (rad)^2);

fprintf('The height of the water in the tank after 10 minutes is %.3f feet \n',height);

fileID = fopen('C:\Users\Owner\Documents\MATLAB\ASEN_2012\Homework_4_Problem_2.txt','w');
fprintf(fileID,'The height of the water in the tank is %.3f feet after a time of ten seconds has ellapsed',height);
fclose(fileID);




