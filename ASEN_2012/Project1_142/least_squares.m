%Purpose: To develop necessary slope and intercept values for a weighted
%least squares line of best fit.
%Input:Time and calorimeter temperature
%Output: m and b values for the input partition of the calorimeter data and the Q matrix associated with those values. 
%Assumptions: Calorimeter specific heat is exact, heat loss from the
%calorimeter is negligible.
%Assigned ID number: 142
%Data Created:10/20/17
%Date(s) modified:10/21/17-10/27/17


function [m, b, Q] = least_squares(time,calor_temp)
%Feed in desired chunk of calor_temp and time interval first in main script
%Calcualte standard deviation of calor_temp in range 1:809
W = zeros(length(time)); %Preallocates weight matrix

A = [time ones(length(time),1)]; %Creates A matrix
d = calor_temp; %Creates d column vector
P = inv(transpose(A)*A)*transpose(A)*d; %calculates a non-weighted least squares fit line
m = P(1); %produces m and b values for line equation y=mx+b
b = P(2);
sigma = uncertainty(time,calor_temp,m,b); %Calculates the uncertainty in values using non-weighted least squares m and b values
for i = 1:length(time) %populates preallocated matrix of zeros with calculated uncertainty 
    W(i,i) = 1/(sigma^2);
end
Q = inv(transpose(A)*W*A); %Calculates Q and new weighted least squares m and b values using created weight matrix
P = inv(transpose(A)*W*A)*transpose(A)*W*d;
m = P(1); %New output m and b values to be used in a weighted least squares fit line y=mx+t+b
b = P(2);



