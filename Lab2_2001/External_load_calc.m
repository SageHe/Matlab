%Read in connectivity and joint coords. Calculate distance between sets of
%joints. Apply mas per length for calculated member length. Apply half of
%calculated force due to member weight at each joint using connectivity of
%joints. Masses are given in Kg

[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput('3PointDesign.txt');
External_loads = zeros(size(joints,1),4);
External_loads(:,1) = 1:size(joints,1);
mag_mass = 1.71e-03; %Kg
joint_mass = .008357; %Kg
mass_per_length = 3.2745e-04; %Kg/in
External_loads(:,4) = mag_mass;
g = 9.81; %m/s^2
for j = 1:size(joints,1)
    num_bars = numel(find(connectivity == j));
    External_loads(j,4) = External_loads(j,4) + (num_bars*joint_mass);
end
for i = 1:size(connectivity,1)
    tmp = connectivity(i,:);
    joint1 = joints(tmp(1),:);
    joint2 = joints(tmp(2),:);
    dist = joint2 - joint1;
    dist = dist.^2;
    dist = sum(dist);
    dist = sqrt(dist);
    mass = dist*mass_per_length;
    mass = mass / 2;
    External_loads([tmp(1) tmp(2)],4) = External_loads([tmp(1) tmp(2)],4) + mass;
end

External_loads(:,4) = External_loads(:,4)*g*-1;


