%Monte Carlo Sim using thermo model varying wind, rocket mass, water mass,
%initial pressure,launch angle, 
%wind +- 2 mph
clear all; close all; clc
N = input('How many tests are you running?\n');
for i = 1:N
    wind_range = [-.8941 .8941];%mph
    wind = 1.7882; %mph
    wind = wind + (wind_range(1) + (wind_range(2) - wind_range(1)).*rand(1,1));
    rock_mass_range = [-.03 .03];
    rock_mass = .122;
    rock_mass = rock_mass + (rock_mass_range(1) + (rock_mass_range(2) - rock_mass_range(1)).*rand(1,1));
    h20_range = [-.02 .02];
    h20_mass = .6;
    h20_mass = h20_mass + (h20_range(1) + (h20_range(2) - h20_range(1)).*rand(1,1));
    pres = 275790; %Pa
    pres_range = [-6894.75 6894.75];
    pres = pres + (pres_range(1) + (pres_range(2) - pres_range(1)).*rand(1,1));
    launch_angle = 40;
    angle_range = [-1 1];
    launch_angle = launch_angle + (angle_range(1) + (angle_range(2) - angle_range(1)).*rand(1,1));
end