%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASEN 2003 - Lab 3: Main Script
%Purpose: 
%   - Imports the data sets into usable variables from LCSDATA.M
%   - Exports the information derived with the data to the LCSMODEL.M
%
%Created: 2/14/2018
%Modified: 2/26/2018
%Creators:
%   - Lucas Zardini
%   - Sage Herrin
%   - Yang Lee
%   - Idam Isnaeni
% --> House - Keeping <--
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables
%radius cm
r = 0.75;
%distance cm
d = 1.53;
%length cm
l = 2.54;
%height cm
h = 2.05;

%% Filename
%ask user about which input file is being used
n = (input('Type the filename you would like to choose: '));

if n == 'T1'
% TEST -> Filename of data
    fname = 'Test1_5pt5V';
    title1 = 'Test 5.5V';
elseif n == 'T2'
% TEST -> Filename of data
    fname = 'Test1_6pt5V';
    title1 = 'Test 6.5V';
elseif n == 'T3'
% TEST -> Filename of data
    fname = 'Test1_7pt5V';
    title1 = 'Test 7.5V';
elseif n == 'T4'
% TEST -> Filename of data
    fname = 'Test1_8pt5V';
    title1 = 'Test 8.5V';
elseif n == 'T5'
% TEST -> Filename of data
    fname = 'Test1_9pt5V';
    title1 = 'Test 9.5V';
elseif n == 'T6'
% TEST -> Filename of data
    fname = 'Test1_10pt5V';
    title1 = 'Test 10.5V';
    
elseif n == 'E1'
% REAL -> Filename of data at 8.5
    fname = 'Test2_8pt5V';
    title1 = '8.5V';
elseif n == 'E2'
% REAL -> Filename of data at 9.5
    fname = 'Test2_9pt5V';
    title1 = '9.5V';
elseif n == 'E3'
% REAL -> Filename of data at 10.5
    fname = 'Test1_10pt5V';
    title1 = '10.5V';
end


%% Input/Output functions
%import LCSDATA
[theta_exp,w_exp,v_exp] = LCSDATA(fname);
%import LCSMODEL
[v_mod] = LCSMODEL(r, d, l, theta_exp, w_exp);

%mm to cm conversion
v_exp = v_exp/100;

%% Error
%difference in model/experimental
error = abs(v_mod - v_exp);
%mean
meanT = mean(error)
%standard deviation
stdT = std(error)
%% Plots
%degrees conversion
theta_exp = theta_exp *(180/pi); %degrees

%%% Velocity vs. angle
%plot MODEL
figure
plot(theta_exp,v_mod)
hold on
%plot EXPERIMENTAL
plot(theta_exp,v_exp)
hold off
xlabel('Theta (degrees)')
ylabel('Velocity (cm/s)')
xlim([theta_exp(1) 6*360])
legend('Model','Experimental')
title('Velocity vs. Angle (' + string(title1) + ')' )

%%% Error plots
figure
plot(theta_exp,error,'r')
title('Error: Velocity vs. Angle (' + string(title1) + ')' )
xlabel('Theta (degrees)')
ylabel('Velocity (cm/s)')
xlim([theta_exp(1) 6*360])