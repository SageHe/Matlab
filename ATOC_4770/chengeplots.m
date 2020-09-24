clear all;close all;clc

data = readtable('Mein0z3.txt');

figure
scatter(data{[1:14],9},data{[1:14],2})
title('Potential temp vs Height 1/3/0Z, Meiningen, Germany')
xlabel('Theta K')
ylabel('Height(m)')

data = readtable('Mein12Z3.txt');

figure
scatter(data{[1:9],9},data{[1:9],2})
title('Potential temp vs Height 1/3/12Z, Meiningen, Germany')
xlabel('Theta K')
ylabel('Height(m)')

data = readtable('Mein0Z4.txt');

figure
scatter(data{[1:18],9},data{[1:18],2})
title('Potential temp vs Height 1/4/0Z, Meiningen, Germany')
xlabel('Theta K')
ylabel('Height(m)')

data = readtable('Mein12Z4.txt');

figure
scatter(data{[1:10],9},data{[1:10],2})
title('Potential temp vs Height 1/4/12Z, Meiningen, Germany')
xlabel('Theta K')
ylabel('Height(m)')



