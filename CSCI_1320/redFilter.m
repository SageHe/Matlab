function [out] = redFilter(current_img, redVal)
b = size(current_img,2);
a = size(current_img,1);
for i = 1:round(b/3) %Counter goes from one to one third the total length of the input image
    for j = 1:a %Counter goes from 1 to the width of the input image
        newimage(j,i,1) = (current_img(j,i,1) * .299) + (current_img(j,i,2) * .587) + (current_img(j,i,3) * .114);
        newimage(j,i,2) = (current_img(j,i,1) * .299) + (current_img(j,i,2) * .587) + (current_img(j,i,3) * .114);
        newimage(j,i,3) = (current_img(j,i,1) * .299) + (current_img(j,i,2) * .587) + (current_img(j,i,3) * .114);
        %Replaces the left third of the image with a grayscale image by
        %multiplying each layer of the new image by the average value of
        %all three layers of the original image
    end
end
for i = round(b/3):round((2*b)/3) %Counter goes from one third the total length of the input image to two thirds the total length of the input image
    for j = 1:a  %Counter goes from 1 to the width of the input image
        newimage(j,i,1) = current_img(j,i,1);
        newimage(j,i,2) = current_img(j,i,2);
        newimage(j,i,3) = current_img(j,i,3);
        %Leaves the middle third of the new image alone, or makes it the
        %exact same as the original input image
    end
end
for i = round(((2*b))/3):b %Counter goes from two thirds the total length of the input image to the full length of the input image
    for j = 1:a %Counter goes from 1 to the width of the input image
        newimage(j,i,1) = (current_img(j,i,1) * redVal) + (current_img(j,i,2) * ((1 - redVal)/2)) + (current_img(j,i,3) * ((1 - redVal)/2));
        newimage(j,i,2) = (current_img(j,i,1) * redVal) + (current_img(j,i,2) * ((1 - redVal)/2)) + (current_img(j,i,3) * ((1 - redVal)/2));
        newimage(j,i,3) = (current_img(j,i,1) * redVal) + (current_img(j,i,2) * ((1 - redVal)/2)) + (current_img(j,i,3) * ((1 - redVal)/2));
        %Applies the red filter to the right third of the image by applying
        %the specified red value input by the user to the red layer of the
        %image, then equally distributing the remaining 'red' left over to
        %the green and blue layers of the image
    end
end
out = newimage;
end
        
        
        
        