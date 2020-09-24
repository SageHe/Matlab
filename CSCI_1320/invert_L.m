function [out] = invert_L(current_img)
for a = 1:size(current_img,1)
    for b = 1:size(current_img,2)
        for c = 1:3
            newImage(a,b,c) =(255 - (current_img(a,b,c)));
        end
    end
end
% This function goes through each column and each row and each layer of the
% matrix and subtracts that value from 255 [max possible matrix element
% value], effectively inverting the colors of the image
out = newImage;
end  