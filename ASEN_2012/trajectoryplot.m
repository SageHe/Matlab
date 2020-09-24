function [out] = trajectoryplot(theta,Trajectory,eccentricity,semimajor_axis)
%This function takes the specified theta vector and eccentricity and plots
%the resulting trajectory of the satellite based on the input
%semimajor-axis and and eccentricity values 
polarplot(theta,Trajectory);
hold on;

legend('Eccentricity = 0, Semimajor-Axis = 1000','Eccentricity = .5, Semimajor-Axis = 1000','Eccentricity = .5, Semimajor-Axis = 1500','Eccentricity = .75, Semimajor-Axis = 1000');
rticklabels('Earth');
end