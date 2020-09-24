%Read in data and format so mean and standard deviation can be found
Data = readtable('HW1data.dat');
Data.Var3 = [];
Data = table2array(Data);
%Calculate mean and standard deviation of each column of data.
m1 = mean(Data(:,1));
std1 = std(Data(:,1));
m2 = mean(Data(:,2));
std2 = std(Data(:,2));