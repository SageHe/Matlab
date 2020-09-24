function [out] = luminance_NL(current_img)
newimage(:,:,1) = ((current_img(:,:,1) * .299) + (current_img(:,:,2) * .587) + (current_img(:,:,3) * .114));
newimage(:,:,2) = ((current_img(:,:,1) * .299) + (current_img(:,:,2) * .587) + (current_img(:,:,3) * .114));
newimage(:,:,3) = ((current_img(:,:,1) * .299) + (current_img(:,:,2) * .587) + (current_img(:,:,3) * .114));
%Sets each layer of the new image equal to the average of all three layers of the
%current image, making a grayscale without using loops
out = newimage;