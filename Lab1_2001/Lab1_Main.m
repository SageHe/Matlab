clear all;close all;clc
File_name = input('What is the name of the text file you want to read? \n');
fileID = fopen(File_name);
Data = textscan(fileID,'%s','Delimiter','\n');
%%Reading and sorting Data from the text file
number_of_forces = str2num(Data{1}{2}(1));
number_of_moments = str2num(Data{1}{2}(4));

force_coords_mat = zeros(number_of_forces,3);
counter1 = 1;
for j = 5:(4 + number_of_forces) %This loop pulls the coordinates of the external forces from the input file
   force_coords_mat(counter1,:) = str2num(Data{1}{j});
   counter1 = counter1 + 1;
end

force_mag_mat = zeros(number_of_forces,4);
counter2 = 1;
for k = (j + 3):(j + 2 + number_of_forces) %This loop pulls the magnitude and direction of the external forces from the input file
    force_mag_mat(counter2,:) = str2num(Data{1}{k});
    counter2 = counter2 + 1;
end

moments_coords = zeros(number_of_moments,3);
counter3 = 1;
for L = (k + 3):(k + 2 + number_of_moments) %This loop pulls the coordinates of the external moments from the input file
    moments_coords(counter3,:) = str2num(Data{1}{L});
    counter3 = counter3 + 1;
end

moments_mag_mat = zeros(number_of_moments,4);
counter4 = 1;
for m = (L + 3):(L + 2 + number_of_moments) %This loop pulls the magnitude and direction of the external moments from the input file
    moments_mag_mat(counter4,:) = str2num(Data{1}{m});
    counter4 = counter4 + 1;
end

support_locations = zeros(6,3);
counter5 = 1;
 for n =(m + 3):(m + 8) %This loops pulls the support coordinates from the input file
    support_locations(counter5,:) = str2num(Data{1}{n});
    counter5 = counter5 + 1;
 end

dummy = zeros((6),4);
force_reactions = zeros((6),4);

for o = (n + 3):(n + 8) %Counts the number of force and moment reactions from the input text file
    if Data{1}{o}(1) == 'F'
        Data{1}{o}(1)='1';
    
    else 
        Data{1}{o}(1)='0';
    end
    dummy(o,:)=str2num(Data{1}{o});
    
end
force_reactions=dummy(n+3:n+8,1:4);
cart_forces = zeros(number_of_forces,3);

for i = 1: number_of_forces %Puts forces into cartesian notation
    cart_forces(i,1) = force_mag_mat(i,1)* (force_mag_mat(i,2)/((force_mag_mat(i,2)^2)+ (force_mag_mat(i,3)^2) + (force_mag_mat(i,4)^2))^(1/2)); %fx
    cart_forces(i,2) = force_mag_mat(i,1)* (force_mag_mat(i,3)/((force_mag_mat(i,2)^2)+ (force_mag_mat(i,3)^2) + (force_mag_mat(i,4)^2))^(1/2)); %fy
    cart_forces(i,3) = force_mag_mat(i,1)* (force_mag_mat(i,4)/((force_mag_mat(i,2)^2)+ (force_mag_mat(i,3)^2) + (force_mag_mat(i,4)^2))^(1/2)); %fz
    

end

moments_by_forces = zeros(number_of_forces,3); %Preallocates number of force matrix of zeros for speed

for i = 1: number_of_forces %Moments of forces
    moments_by_forces(i,:) = cross(force_coords_mat(i,:),cart_forces(i,:));
end

couple_moments = zeros(number_of_moments,3); %Preallocates a couple moments vector of zeros for speed

for i = 1:number_of_moments %Getting the magnitude into the x, y, and z directions
    couple_moments(i,1) = moments_mag_mat(i,1) *(moments_mag_mat(i,2)/((moments_mag_mat(i,2)^2)+(moments_mag_mat(i,3)^2)+(moments_mag_mat(i,4)^2))^(1/2));
    couple_moments(i,2) =moments_mag_mat(i,1) *(moments_mag_mat(i,3)/((moments_mag_mat(i,2)^2)+(moments_mag_mat(i,3)^2)+(moments_mag_mat(i,4)^2))^(1/2));
    couple_moments(i,3) = moments_mag_mat(i,1) *(moments_mag_mat(i,4)/((moments_mag_mat(i,2)^2)+(moments_mag_mat(i,3)^2)+(moments_mag_mat(i,4)^2))^(1/2));
end

%% Solving matrix for unknown magnitudes
sum_of_forces = sum(cart_forces,1); %Sums the force vectors
if number_of_moments ~= 0
    sum_of_moments = couple_moments + sum(moments_by_forces); %Sums all the moments
else 
    sum_of_moments = sum(moments_by_forces,1);
end

v1 = cat(1,sum_of_forces,sum_of_moments);
v1 = horzcat(sum_of_forces,sum_of_moments);
v1 = rot90(v1,3); %Puts the sum of the forces in the x,y, and z directions in a vector
force_counter = 0;
moment_counter = 0;

unit_reactions = force_reactions;

for i = 1:length(force_reactions) %Makes the reaction forces unit vectors
    
    unit_reactions(i,2) = (force_reactions(i,2)) /( (force_reactions(i,2)^2) + (force_reactions(i,3)^2)+(force_reactions(i,4)^2))^(1/2);
    unit_reactions(i,3) = (force_reactions(i,3)) /( (force_reactions(i,2)^2) + (force_reactions(i,3)^2)+(force_reactions(i,4)^2))^(1/2);
    unit_reactions(i,4) = (force_reactions(i,4)) /( (force_reactions(i,2)^2) + (force_reactions(i,3)^2)+(force_reactions(i,4)^2))^(1/2);
   
end

matrix = zeros(6,6);

for i = 1:length(matrix) %This loop creates the 6x6 A matrix
    if force_reactions(i,1) == 1
        matrix(1,i) = unit_reactions(i,2); 
        matrix(2,i) = unit_reactions(i,3);
        matrix(3,i) = unit_reactions(i,4);
        
        moments_by_supports = cross(support_locations(i,:),unit_reactions(i,2:end));
        
        matrix(4,i) = moments_by_supports(1);
        matrix(5,i) = moments_by_supports(2);
        matrix(6,i) = moments_by_supports(3);
    else
        matrix(4,i) = unit_reactions(i,2);
        matrix(5,i) = unit_reactions(i,3);
        matrix(6,i) = unit_reactions(i,4);
        
    end
    
end

v1 = v1 .* (-1); 

X = inv(matrix)*v1; %outputs the X matrix that has been solved for 
disp(X);

dlmwrite('Example_output',X,'delimiter','\n')
%Got help from TAs and peers  