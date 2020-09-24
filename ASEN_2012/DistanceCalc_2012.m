clear all; close all;
X1 = input('Enter the x-coordinate of the first point:');
Y1 = input('Enter the y-coordinate of the first point:');
X2 = input('Enter the x-coordinate of the second poin:');
Y2 = input('Enter the y-coordinate of the second point:');

distance = sqrt((X2 - X1).^2 + (Y2 - Y1).^2);

fprintf('The distance between the two points is %f units \n',distance);