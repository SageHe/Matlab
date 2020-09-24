clear all;close all;clc

t1 = 10;
t2 = 20;
t3 = 30;
% syms t1 t2 t3

DCM321 = [cosd(t2)*cosd(t1) cosd(t2)*sind(t1) -sind(t2); sind(t3)*sind(t2)*cosd(t1)-cosd(t3)*sind(t1) sind(t3)*sind(t2)*sind(t1)+cosd(t3)*cosd(t1) sind(t3)*cosd(t2);cosd(t3)*sind(t2)*cosd(t1)+sind(t3)*sind(t1) cosd(t3)*sind(t2)*sind(t1)-sind(t3)*cosd(t1) cosd(t3)*cosd(t2)]

DCM313 = [cosd(t3)*cosd(t1)-sind(t3)*cosd(t2)*sind(t1) cosd(t3)*sind(t1)+sind(t3)*cosd(t2)*cosd(t1) sind(t3)*sind(t2); ...
            -sind(t3)*cosd(t1)-cosd(t3)*cosd(t2)*sind(t1) -sind(t3)*sind(t1)+cosd(t3)*cosd(t2)*cosd(t1) cosd(t3)*sind(t2); ...
            sind(t2)*sind(t1) -sind(t2)*cosd(t1) cosd(t2)]
