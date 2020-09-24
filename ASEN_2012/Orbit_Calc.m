function [theta,r] = Orbit_Calc(semimajor_axis,eccentricity)
%This funtion takes the input a and e values given by the user and
%implements them in the equation to calcualte the trajectory of the
%satellite
theta = 0:.01:2*pi;

r = semimajor_axis./(1 - eccentricity.*cos(theta));


