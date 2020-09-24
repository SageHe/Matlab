function [out] = Famous_Me(current_img, location)
newimage = location; %preallocates the current image to be equal to the location image
for i = 1:size(current_img,1) %counter that goes from one to the length of the current image, which in this case is the small version of the picture of me
    for j = 1:size(current_img,2) %Counter that goes from one to the width of the current image, which in this case is the small version of the picture of me
            if current_img(i,j,:) >= 1 %If the value of the current image in any layer is andything but absolutely black, replace that pixel in the location image with with appropriate image from the picture of myself
                newimage(i + 100,j + 100,:) = current_img(i,j,:);
            end
    end
end
out = newimage;
end
            