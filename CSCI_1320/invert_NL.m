function [out] = invert_NL(current_img)
newimage(:,:,1) = 255 - (current_img(:,:,1));%Adds the value of the entered brightess to the entire first level of the image
newimage(:,:,2) = 255 - (current_img(:,:,2));%Adds the value of the entered brightess to the entire second level of the image
newimage(:,:,3) = 255 - (current_img(:,:,3));%Adds the value of the entered brightess to the entire third level of the image
out = newimage;%Makes the new levels that have been altered and and combines them to make a new image 
