%% Title
%Code Written By: Adam Boylston and Sage Herrin 
%Date Created: 10/30/17
%Last Modified: 10/30/2017

%% Housekeeping
close all
clear
clc

%% Problem 5.1.1
%Reads boulder image and converts to double
boulder = double(imread('boulder.jpg'));

boulder_bright = boulder + 100;

%Display new image side by side with old image
subplot(1,2,1),imshow(boulder,[0,255])
title('Original')
subplot(1,2,2),imshow(boulder_bright,[0,255])
title('Increased Exposure')

pause
close

%% Problem 5.1.2
%creates arbitrary A matrix (4x4 grayscale)
A = [50 100 150 200;50 100 150 200;50 100 150 200;50 100 150 200];
E = [0 0 0 1; 0 1 0 0; 0 0 1 0; 1 0 0 0];

subplot(1,3,1),imshow(A,[0,255])
title('A')
subplot(1,3,2),imshow(A*E,[0,255])
title('A*E')
subplot(1,3,3),imshow(E*A,[0,255])
title('E*A')

pause
close

%% Problem 5.1.3
%Reads photo1 and converts to double
photo1 = double(imread('photo1.jpg'));

%sets variable equal to size of image for indexing
m = length(photo1);
hshift_mag = 140;

photo1_hshift = zeros(m);
photo1_eye = eye(m);
photo1_hshift(:,1:hshift_mag) = photo1_eye(:,(m-(hshift_mag-1):m));
photo1_hshift(:,hshift_mag+1:m) = photo1_eye(:,(1:(m-hshift_mag)));
photo1_shifted = photo1*photo1_hshift;

subplot(1,3,1),imshow(photo1,[0,255])
title('Original')
subplot(1,3,2),imshow(photo1_shifted,[0,255])
title('Shifted Horizontally 140px')
subplot(1,3,3),spy(photo1_hshift)
title('Shift Matrix Spy')

pause
close

%% Problem 5.1.4
%Performs vertical and horizontal shifts on the input image.
vshift_mag = 100;
photo1_hvshift = zeros(m);

photo1_hvshift(1:vshift_mag,:) = photo1_eye(m-(vshift_mag-1):m,:);
photo1_hvshift(vshift_mag+1:m,:) = photo1_eye(1:m-vshift_mag,:);

subplot(1,3,1),imshow(photo1,[0,255])
title('Original')
subplot(1,3,2),imshow(photo1_hvshift*photo1_shifted,[0,255]) %switch order multiplied
title('Shifted Horizontally 140px, Vertically 100px')
subplot(1,3,3),spy(photo1_hvshift*photo1_hshift) 
title('Shift Matrix Spy')

pause
close

%% Problem 5.1.5
%Determine what matrix is required to flip an image upside down.

%Create a diagonal matrix equal in size to photo1

updown = eye(length(photo1));

%Make B a reversed diagonal matrix

updown = updown([256:-1:1],:);

%Multiply B by the the photo 'photo1', paying close attention to the order
%of multiplication 

temp = updown * photo1; %Temporary variable of manipulated photo1
imshow(temp,[0 255])

subplot(1,3,1),imshow(photo1,[0 255])
title('Original')
subplot(1,3,2),imshow(temp,[0 255])
title('Inverted')
subplot(1,3,3),spy(updown)
title('Identity matrix to invert image')

pause
close

%% Problem 5.2.1
%Create a DCT (Discrete Cosine Transform) matrix and verify that for 5X5
%DCT matrix called C, C = C^-1, or show that C^2 = I5 (5X5 identity matrix)
mat = ones(5);
DCT_mat = DCT(mat);

%Verify that DCT_mat = DCT_mat^-1 => DCT_mat*DCT_mat = I5(5X5) identity
%matrix by multiplying DCT_mat by itself

temp = DCT_mat * DCT_mat;

disp(temp)

%% Problem 5.2.2
%Compute determinant of DCT matrix for 3 value of n. 

mat = ones(3);

DCT_mat = DCT(mat);

%Check determinant of DCT_matrix and speculate on relationship and what it
%says about C^-1 

check = det(DCT_mat);

disp(check)

%Determinant changes sign after 2 increases of n, e.g. determinant of 3 and
%4 = -1, determinant of 5 and 6 = 1, determinant o 7 and 8 = -1, etc. 

%% Problem 5.2.3
%Use DCT matrix code to to compute eigenvalues of C for n = 512, use MATLAB
%function 'eig'. Plot eigenvalues on complex plane. Speculate on pattern in
%the determinants given knowledge of eigenvalues. 

DCT_mat = DCT(ones(512)); %512X512 matrix of ones

temp = eig(DCT_mat); %Computes eigen values of DCT matrix

plot(temp,'*') %Plot Eigenvalues, represented by asterisks.


%% Problem 5.2.4
DCT = zeros(size(photo1,1));
for i = 1:size(photo1,1) %DCT matrix build
    for j = 1:size(photo1,2)
        DCT(i,j) = sqrt((2/length(photo1)))*cos((pi *(i - (1/2))*(j - (1/2)))/length(photo1));
    end
end

%given two dimensional transform
two_d = DCT * photo1 * DCT;

%showing that C^-1 = C, since multiplying the 2D transform by either
%restores the original image
restored1 = DCT*two_d*DCT;
restored2 = inv(DCT)* two_d * inv(DCT);

subplot(1,3,1),imshow(two_d,[0 255])
title('2-D')
subplot(1,3,2),imshow(restored1,[0 255])
title('Restored')
subplot(1,3,3),imshow(restored2,[0 255])
title('Restored')

pause
close

%% Problem 5.2.5
%Uses the DCTfunc to compress images using various values of p
subplot(1,4,1),imshow(DCTfunc(boulder,1),[0 255])
title('Original, p = 1')
subplot(1,4,2),imshow(DCTfunc(boulder,.7),[0 255])
title('Restored, p = 0.7')
subplot(1,4,3),imshow(DCTfunc(boulder,.3),[0 255])
title('Restored, p = 0.3')
subplot(1,4,4),imshow(DCTfunc(boulder,.05),[0 255])
title('Restored, p = 0.05')

pause
close

%% Problem 5.2.7
%Determine compression ratio for p = .5,.3, and .1. Use number of nonzero
%entries in transformed image matrix Y as sub for file size. What value of
%p gives best CR while maintaining reasonable image quality.

[X,Q1] = DCTfunc(boulder,.5);

[Y,Q2] = DCTfunc(boulder,.3);

[Z,Q3] = DCTfunc(boulder,.1);

%Use compressionratio as number of nonzero elements in uncompressed image
%over number of nonzero elements in compressed image. 

CR_5 = (nnz(boulder) / nnz(Q1));
CR_3 = (nnz(boulder) / nnz(Q2));
CR_1 = (nnz(boulder) / nnz(Q3));

fprintf('Compression ratio for a p value of .5 is %.4f \n', CR_5);
fprintf('Compression ratio for a p value of .3 is %.4f \n', CR_3);
fprintf('Compression ratio for a p value of .1 is %.4f \n', CR_1);
%Looks like p value of approximately .15 yields best compression ratio
%while maintaining good image quality.


%% Problem 5.2.8
%Read in boulder_noisy.jpg and compute Y matrix for several values of p.
%Speculate on how much DCT helps.

boulder_noisy = double(imread('boulder_noisy.jpg'));
%Use the DCTfunc to reduce noise of boulder_noisy uisng various p values
dct_1 = DCTfunc(boulder_noisy,.7);
dct_2 = DCTfunc(boulder_noisy,.5);
dct_3 = DCTfunc(boulder_noisy,.3);

subplot(1,4,1),imshow(boulder_noisy, [0 255])
title('Original')
subplot(1,4,2),imshow(dct_1, [0 255])
title('DCT .7 P value')
subplot(1,4,3),imshow(dct_2, [0 255])
title('DCT .5 P value')
subplot(1,4,4),imshow(dct_3, [0 255])
title('DCT .3 P value')

pause
close
%DCT function help to clean up noise a little bit, but only to a certain
%extent. Any P value below approximately .55 the image starts to become to
%compressed. P value of approximately .57 yields optimal noise reduction.
DCTslider(boulder)
function DCTslider(boulder)
handles.Image = boulder;
handles.fig = figure;
handles.axes1 = axes('Units','pixels','Position',[50 100 400 400]);
handles.slider = uicontrol('Style','slider','Position',[50 50 400 20],'Min',0,'Max',1);
set(handles.slider,'Value',1);
handles.Listener = addlistener(handles.slider,'Value','PostSet',@(s,e) DCTplot);
imshow(handles.Image,'Parent',handles.axes1);
guidata(handles.fig, handles);
imshow(DCTfunc(boulder,1),[0 255])
bl1 = uicontrol('Style','text','Position',[45,25,25,25],...
                'String','0');
bl2 = uicontrol('Style','text','Position',[430,25,25,25],...
                'String','1');
bl3 = uicontrol('Style','text','Position',[200,25,100,25],...
                'String','P = 1');

function DCTplot
    slider_value = get(handles.slider,'Value');
    h = DCTfunc(boulder,slider_value);
    %handles.Image=imfilter(handles.Image,h,'conv');
    axes(handles.axes1);
    imshow(h,[0 255])
    caption = sprintf('P = %.2f', slider_value);
    bl3 = uicontrol('Style','text','Position',[200,25,100,25],...
                'String',caption);
end
end

function [boulder_restored,Q] = DCTfunc(boulder,p)
n = size(boulder,1);
dct = DCT(ones(length(boulder)));
temp = dct * boulder * dct;
for i = 1:n
    for j = 1:n
        if i + j > p*2*n
            temp(i,j) = 0;
%         else
%             C(i,j) = sqrt((2/length(boulder)))*cos((pi *(i - (1/2))*(j - (1/2)))/length(boulder));
        end
    end
end
%boulder_2D = C * boulder * C;
Q = temp;
boulder_restored = dct * temp * dct;
end
