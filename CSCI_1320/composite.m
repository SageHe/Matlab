function [out] = composite (current_img)
orig_img = current_img;
red_img(:,:,1) = current_img(:,:,1);
red_img(:,:,2) = zeros(size(current_img,1),size(current_img,2));
red_img(:,:,3) = zeros(size(current_img,1),size(current_img,2));
%Makes a red philter of the image by making the red layer of 'red_img' equal to
%the origianl red layer of the image and making the green and blue layers
%of the 'red_img' equal to a zeros matrix of equal size to the original
%image
green_img(:,:,2) = current_img(:,:,2);
green_img(:,:,1) = zeros(size(current_img,1),size(current_img,2));
green_img(:,:,3) = zeros(size(current_img,1),size(current_img,2));
%Makes a green philter of the image by making the green layer of 'green_img' equal to
%the origianl green layer of the image and making the red and blue layers
%of the 'green_img' equal to a zeros matrix of equal size to the original
%image
blue_img(:,:,3) = current_img(:,:,3);
blue_img(:,:,1) = zeros(size(current_img,1),size(current_img,2));
blue_img(:,:,2) = zeros(size(current_img,1),size(current_img,2));
%Makes a blue philter of the image by making the blue layer of 'blue_img' equal to
%the origianl blue layer of the image and making the red and green layers
%of the 'blue_img' equal to a zeros matrix of equal size to the original
%image
a = size(current_img,1);
b = size(current_img,2);
for i = 1:a
    for j = 1:b
        new_img(i,j,:) = current_img(i,j,:);
        new_img(i + a,j,:) = red_img(i,j,:);
        new_img(i,j + b,:) = green_img(i,j,:);
        new_img(i + a,j + b,:) = blue_img(i,j,:);
    end
end
%Goes through using a counter that goes from 1 to the length and width of
%the original image and, based on the quadrant of the new image specified
%by adding the counter to the length or width of the original image, makes
%each corner of the new image either the original image, the red layer of
%the image, the green layer of the image, or the blue layer of the image
out = new_img;
end

