clear all;close all;clc

data = xlsread('Testdata.xlsx');

test_no = data(:,1);
F = data(:,2); % Force
a = data(:,3); % distance from left most support 
w = data(:,4); % width of beam top down view
d_f = data(:,5); % distance from the failure location of beam to closet support
