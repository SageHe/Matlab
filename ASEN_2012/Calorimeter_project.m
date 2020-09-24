%Purpose: To determine the specific heat and its uncertainty of an unknown substance using best
%fit line methods and error propagation and analysis techniques.
%Input:Time, calorimeter temperature, and water temperature in one Excel
%file, calorimeter and sample mass and their uncertainties, calorimeter
%specific heat.
%Output: Sample specific heat and its uncertainty, t0,t1, and t2 values and
%their associated uncertainties. 
%Assumptions: Calorimeter specific heat is exact, heat loss from the
%calorimeter is negligible.
%Assigned ID number: 142
%Data Created:10/20/17
%Date(s) modified:10/21/17-10/27/17

clear all;close all;clc
%For excel Data, first column is time, second is calorimeter temp, third is
%water temp. Initial water temp is same as initial sample temp.

Data = xlsread('CalorimeterData.xlsx');

time = Data(:,1); %First column of Excel data, time
calor_temp = Data(:,2); %Second column of Excel data, calorimeter temperature
water_temp = Data(:,3); %Third column of Excel data, water temperature

calor_mass = 313.5; %Mass of the calorimeter in grams
calor_mass_uncert = .05; %Uncertinty in Mass of calorimeter in grams
calor_C = .214; %Specific heat of the calorimeter in calories/gram*C, assumed exact
sample_mass = 91.75; %Mass of sample in grams
mass_uncert = .05; %Uncertainty in the mass of the sample in grams

%By inspection of the graph it appears that the temp rapidly starts to
%increase at the 809th data point => t0 = calor_temp(809)

%Fit the first 809 data points with a weighted least squares line
[m1, b1] = least_squares(time(1:809),(calor_temp(1:809)));
%Calculates uncertainty in line fit to specified portion of data (see
%comment in function)
sigmat0 = uncertainty(time(1:809),(calor_temp(1:809)),m1,b1);
x1 = time(1:809); %First best fit line
y1 = m1*x1 + b1;
%Calculates the initial temperature of the calorimeter 
t0 = (m1*time(809) + b1); 
%Fit data from the 1057th data point to the last data point with a weighted
%least squares line
[m2, b2, Q2] = least_squares(time(1057:end),(calor_temp(1057:end)));
x2 = time(809:end); %Last best fit line
x3 = time(1057:end); %Last best fit line
y2 = m2*x2 + b2;
y3 = m2*x3 + b2;
t_low = t0; %Assignes t_low as t0 since they are the same, specified in assignment
t_high = y2(1); %Point on extrapolated line corresponding to same t value as t0

t_avg = (t_high + t_low) / 2; %Finds average t on middle portion of graph where the sample is added to the calorimeter
%Fit data from point 840 to 922 with weighted least squares line 
[m3, b3] = least_squares(time(810:922),(calor_temp(810:922)));
x4 = time(810:922);%Middle best fit line
y4 = m3*x4+b3;

t = ((t_avg - b3)/m3); %solves for time value t in order to solve for t2

t1 = mean(water_temp(1:809)); %Average water temperature in first 809 data points
sigmat1 = std(water_temp(1:809)); %Standard deviatio in those data points

t2 = m2*t + b2; %Calculates t2
sigmat2 = sqrt([t 1]*Q2*[t;1]); %Calculates uncertianty in t2 using principles in calculating incertainty in an extrapolated point of data from lecture notes, lecture 12
%calculate value of Cs given equation one
Cs = calor_mass*calor_C*(t2 - t0)/(sample_mass*(t1 - t2));

%Convert calorimeter specific heat to J/g*C
Cs = Cs * 4.184; %Using the conversion of 4.184 joules per calories

%Find the uncertainty in the final Cs calculation using partial
%derivatives, given uncertainties, and calculated uncertainties

partialWRS_Mc = (calor_C*(t2 - t0)) / (sample_mass*(t1 - t2));
partialWRS_Ms = -(calor_mass*calor_C*(t2 - t0)) / ((sample_mass)^2*(t1 - t2));
partialWRS_t0 = -(calor_mass*calor_C) / (sample_mass*(t1 - t2));
partialWRS_t1 = -(calor_mass*calor_C*(t2 - t0)) / ((t1 - t2)^2*sample_mass);
partialWRS_t2 = (calor_mass*calor_C*(t1 - t0)) / (sample_mass*(t2 - t1)^2);

sigma_Cs = abs(partialWRS_Mc)*calor_mass_uncert + abs(partialWRS_Ms)*mass_uncert + abs(partialWRS_t0)*sigmat0 + abs(partialWRS_t1)*sigmat1 + abs(partialWRS_t2)*sigmat2;
sigma_Cs = sigma_Cs * 4.184; %Conversion from calories to joules in the uncertainty
hold on
plot(Data(:,1),Data(:,2)) %Plots the temperature of the calorimeter vs time
plot(x1,y1) %Least squares-fit line for first third of data
plot(x2,y2,'--') %Leas squares-fit line extrapolated backwards to find t2
plot(x3,y3) %Least squares-fit line of third portion of data
plot(x4,y4) %Least squares fit of middle portio of data
title('Calorimeter Temperature vs Time')
xlabel('Time(s)')
ylabel('Temperature(Celsius)')
legend('Calorimter Temperature','Least Squares-Fit for calorimeter initial temp','Extended Least Squares-fit of sample final temp for extrapolation','Least Squares-fit line for sample final temp','Least Squares-fit line for calorimeter-sample temp before equilibrium')
hold off
figure 
plot(Data(:,1),Data(:,3))
title('Water temperature')
xlabel('Time(s)')
ylabel('Temperature(Celsius)')
legend('Water Temperature')

