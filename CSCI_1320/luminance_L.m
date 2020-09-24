function [out] = luminance_L(current_img)
for a = 1:size(current_img,1) %Sets a counter a from 1 to the size of the rows of the current image
    for b = 1:size(current_img,2) %Set a counter b from 1 to the size of the columns of the current image
        newimage(a,b,1) = ((current_img(a,b,1) * .299) + (current_img(a,b,2) * .587) + (current_img(a,b,3)));
        newimage(a,b,2) = ((current_img(a,b,1) * .299) + (current_img(a,b,2) * .587) + (current_img(a,b,3)));
        newimage(a,b,3) = ((current_img(a,b,1) * .299) + (current_img(a,b,2) * .587) + (current_img(a,b,3)));
        %Makes every layer of the new image equal to the average of all
        %three of the layers of the original image, making a grayscale
        %image
    end
end
out = newimage;