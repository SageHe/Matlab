function [out] = flag(current_img1, current_img2, current_img3)
b = size(current_img1,1);
a = size(current_img1,2);
for x = 1:size(current_img2,1) %Counter goes from one to the length of the image
    for y = 1:size(current_img2,2) %Counter goes from on the the width of the image
        for z = 1:3 %Counter goes from one to three
            modcurrent_img2(x,y,z) = (255 - (current_img2(x,y,z)));
            %Creates an inverted version of the orginal image saved as
            %'current_img2'
        end
    end
end
modcurrent_img1(:,:,1) = current_img1(:,:,1);
modcurrent_img1(:,:,2) = 0;
modcurrent_img1(:,:,3) = 0;
%Makes the image saved as 'current_img1' a red philter version of that
%image
modcurrent_img3(:,:,3) = current_img3(:,:,3);
modcurrent_img3(:,:,1) = 0;
modcurrent_img3(:,:,2) = 0;
%Makes the image saved as 'current_img2' a blue philter version of that
%image
for i = 1:a
    for j = 1:b
        newimage(i,j,:) = modcurrent_img1(i,j,:);
        newimage(i + b,j,:) = modcurrent_img2(i,j,:);
        newimage(b + b + i,j,:) = modcurrent_img3(i,j,:);
        %Puts all three of the modified images together
    end
end
out = newimage;
end
        
        