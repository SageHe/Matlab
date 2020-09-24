function [out] = binaryMask(current_img)
newimage = current_img; %Preallocates the newimage to be equal in size to the current image
for i = 1:size(current_img,1) %Counter goes from one to the length of the input image
    for j = 1:size(current_img,2) %Counter goes from one to the width of the input image
        if current_img(i,j,:) <= 100 %If the value of any layer is less than or equal to 100, replace that value with 255, or make it white
            newimage(i,j,:) = 255;
        elseif current_img(i,j,:) >=100 %If the value of any layer is greater than or equal to 100, replace that value with 0, or make it black
            newimage(i,j,:) = 0;
        end
    end
end
out = newimage;
end


