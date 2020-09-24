function [out] = addRandomNoise_NL(current_img)
a = size(current_img);%Creates a matrix of the size of the image
b = randi(([-255 255]),a);%For every element of a it creates a random integer between -255 and 255
c = uint8(b);%Type casts the b matrix to uint8
out = current_img + c;%Adds the value of the elements of the matrix of c to the current image matrix
end