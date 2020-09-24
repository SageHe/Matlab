%Purpose: To calculate the uncertainty in measurements
%Input:Time,calorimeter temperature, m, and b values for weight line of
%best fit for specified partition of data
%Output: Uncertainty sigma of input values
%Assumptions: N/A
%Assigned ID number: 142
%Data Created:10/20/17
%Date(s) modified:10/21/17-10/27/17


function out = uncertainty(time,calor_temp,m,b)
sigma = [];
for i = 1:length(calor_temp)
    sigma = [sigma (calor_temp(i) - b - m*time(i)).^2]; 
end
sigma = sum(sigma);
sigma = sqrt(1/(length(time) - 2)*sigma);  %Calculates the uncertainty in measured value using equation 8.15 from John R. Taylor's An Introduction to Error Analysis 
out = sigma;