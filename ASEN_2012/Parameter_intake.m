function [semimajor_axis,eccentricity] = Parameter_intake()
%This function prompts the user to enter the two required parameters for
%the calcualation of the orbit of the spacecraft
semimajor_axis = input('Enter the value of the semimajor axis of the ellipse that represents the orbit of the satellite \n');
eccentricity = input('Enter the value of the eccentricity of the ellipse that represents the orbit of the satellite \n');
while (eccentricity <0) || (eccentricity >= 1)
    eccentricity = input('Please enter an eccentricity value greater than zero and less than one \n');
end
end

