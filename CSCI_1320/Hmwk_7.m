% This script creates an image processing menu driven application

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% your information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;

% Display a menu and get a choice
choice = menu('Choose an option', 'Exit Program', 'Load Image', ...
    'Display Image', 'Invert Image With Loops', 'Invert Image No Loops', ...
    'Add Random Noise With Loops', 'Add Random Noise No Loops', 'Luminance With Loops', ...
    'Luminance No Loops', 'composite', 'color swap', 'red filter', 'flag', 'binary mask', 'famous me');  % as you develop functions, add buttons for them here
 
% Choice 1 is to exit the program
while choice ~= 1
   switch choice
       case 0
           disp('Error - please choose one of the options.')
           % Display a menu and get a choice
           choice = menu('Choose an option', 'Exit Program', 'Load Image', ...
    'Display Image', 'Brighten Image');  % as you develop functions, add buttons for them here
        case 2
           % Load an image
           image_choice = menu('Choose an image', 'lena1', 'lena1_small', 'wrench', 'mandril1', 'yoda', 'sully', 'shrek', 'yoda_small', 'redBaloon');
           switch image_choice
               case 1
                   filename = 'lena1.jpg';
               case 2 
                   filename = 'lena1_small.jpg';
               case 3
                   filename = 'wrench1.jpg'; 
               case 4
                   filename = 'mandrill1.jpg';
               case 5
                   filename = 'yoda.bmp';
               case 6
                   filename = 'sully.bmp';
               case 7 
                   filename = 'shrek.bmp';
               case 8
                   filename = 'yoda_small.bmp';
               case 9
                   filename = 'redBaloon.jpg';
               % fill in cases for all the images you plan to use
           end
           current_img = imread(filename);
       case 3
           % Display image
           figure
           imagesc(current_img);
           if size(current_img,3) == 1
               colormap gray
           end
           
       case 4
           % Invert Image using loops
           % 2. Call the appropriate function
             newImage = invert_L(current_img);
           % 3. Display the old and the new image using subplot
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           % Displays the original image and the new image side by side
           % using subplot
           % 4. Save the new image to a file
           imwrite(newImage,'Alteredimage1.jpg')
           
              
       case 5
           %Invert Image No loops
           %Call  the function to invert the image without using loops
           newImage = invert_NL(current_img);
           %Display the old and the new image using subplot
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image to a file
           imwrite(newImage,'Alteredimage2.jpg')
       case 6
           %Add random noise to the image while using loops
           %Call the function to add random noise using loops
           newImage = addRandomNoise_L(current_img);
           %Display original image and new image together using subplot
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image to a file
           imwrite(newImage,'Alteredimage3.jpg')
       case 7
           %Add random noise to the image without using loops
           %Call the function to add random noise without using loops
           newImage = addRandomNoise_NL(current_img);
           %Display original image and new image together using sublot
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image to a file
           imwrite(newImage,'Alteredimage4.jpg')
       case 8
           %Make the color image a grayscale image while using loops
           %Call the function to make the image grayscale using loops
           newImage = luminance_L(current_img);
           %Display the current image and the new image together using
           %subplot
           subplot(1,2,1), imshow(current_img)
           subplot(1,2,2), imshow(newImage)
           %Save the new new image to a file
           imwrite(newImage,'Alteredimage5.jpg')
       case 9
           %Make the color image a grayscale image without using loops
           %Call the function to make the image grayscale without using
           %loops
           newImage = luminance_NL(current_img);
           %Display the current image and the new image together using
           %subplots
           subplot(1,2,1), imshow(current_img)
           subplot(1,2,2), imshow(newImage)
           %Save the new image to a file
           imwrite(newImage,'Alteredimage6.jpg')
       case 10
           %Make a composite image consisting of the original image in the
           %top left corner and a red, green, and blue image in the other
           %three corners
           %Call the function to make the composite image described above
           newImage = composite(current_img);
           %Display orginal image and new image side by side using subplots
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image to a file
           imwrite(newImage,'Alteredimage7.png')
       case 11
           %Replace pixels of one color with pixels of another color with some leniency
           %Call the function to swap pixels of one color with another
           %color
           r1 = input('please input the target red value of the original image\n');
           g1 = input('please input the target green value of the original image\n');
           b1 = input('please input the target blue value of the original image\n');
           allowed = input('Please enter the desired value lenience\n');
           r2 = input('Please input the new red value of the new image\n');
           g2 = input('Please input the new green value of the new image\n');
           b2 = input('Please input the new blue value of the new image\n');
           newImage = color_swap(current_img,r1,g1,b1,allowed,r2,g2,b2);
           %Display original image and new image side by side using subplot
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image as a file
           imwrite(newImage,'Alteredimage8.png')
       case 12 
           %Apply a red filter to the desired image
           %Call the function that will apply the red filter to the right
           %third of the picture, make the left third a grayscale image,
           %and leaves the middle alone
           redVal = input('please input the red filter value between 0 and 1\n');
           newImage = redFilter(current_img, redVal);
           %Display the original image and the new image side by side using
           %subplots
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the new image as a file
           imwrite(newImage,'Alteredimage9.png')
       case 13
           %Make a flag collage out of at least three different images
           %Call the function to make the flag collage
           current_img1 = imread('Elrond.jpeg');
           current_img2 = imread('Sauron_eye.jpg');
           current_img3 = imread('OG_Thranduil.png');
           newImage = flag(current_img1, current_img2, current_img3);
           %Display the resulting flag image
           imagesc(newImage)
           %Save the flag image to a file
           imwrite(newImage,'Alteredimage10.png')
       case 14
           %Make a binary mask of an image
           %Call the function that makes the binary mask of the image
           newImage = binaryMask(current_img);
           %Display the two images side by side using subplots
           subplot(1,2,1), imagesc(current_img)
           subplot(1,2,2), imagesc(newImage)
           %Save the image to a file
           imwrite(newImage,'Alteredimage11.png')
       case 15
           %Takes a picture of a famous area, and puts a smaller version of
           %yourself onto that image
           %Call the function to put yourself on to the famous location
           current_img = imread('sage_dark.jpg');
           small_img = imresize(current_img,.25);
           location = imread('The_Shire.jpg');
           newImage = Famous_Me(small_img, location);
           %Display the background, you, and the final image side by side
           %using subplots
           subplot(1,3,1), imagesc(location)
           subplot(1,3,2), imagesc(small_img)
           subplot(1,3,3), imagesc(newImage)
           %Save the final image to a file
           imwrite(newImage,'Alteredimage12.png')
   end
   % Display menu again and get user's choice
   choice = menu('Choose an option', 'Exit Program', 'Load Image', ...
    'Display Image', 'Invert Image With Loops', 'Invert Image No Loops', ...
    'Add Random Noise With Loops', 'Add Random Noise No Loops', 'Luminance With Loops', ...
    'Luminance No loops', 'composite', 'color swap', 'red filter', 'flag', 'binary mask', 'famous me');  % as you develop functions, add buttons for them here
end
