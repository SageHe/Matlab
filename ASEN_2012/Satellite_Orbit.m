%This script uses three sub-functions to take in the a and e values the
%user wants to use to calculate the trajectory of the satellite, calculate
%the trajectory based on the input values, and plot the calculated
%trajectory

%Call the function that prompts the user to input the desired a and e
%values

[semimajor_axis,eccentricity] = Parameter_intake;

%Call the function that takes in the a and e values given by the user and
%implements them into the equation to calculate the orbit of the satellite

[theta,Trajectory] = Orbit_Calc(semimajor_axis,eccentricity);

%Call the function that takes the output trajectory of the satellite and
%plots it with Earth at the center

trajectoryplot(theta,Trajectory,eccentricity,semimajor_axis)

