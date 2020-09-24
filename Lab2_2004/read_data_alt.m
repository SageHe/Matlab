% function [Isp, delta_v, data] = read_data_alt()
clear all;close all;clc
folder_name = input('What folder are you reading out of?\n');
D = dir(folder_name);
Isp = [];
delta_v = [];
fuel_mass = .8;
for i = 3:size(D,1)
data = xlsread(D(i).name);
new_data = data(:,3);
new_data = convforce([new_data],'lbf','N');
time = linspace(0,length(new_data)/1652,12000);
time = time';
plot(time,new_data)
title('Thrust Vs Time(600g)')
xlabel('Time (s)')
ylabel('Thrust (N)')
% A = diff(new_data);
ind1 = find(new_data > 10);
ind1 = ind1(1);
% ind1 = ind1/1652;
trunc = find(new_data == max(new_data));
temp = new_data(trunc:end);
min_range = find((new_data > -28) & (new_data < -24));
ind2 = (min_range(round(numel(min_range)/2)));
% ind2 = ind2/1652;
% ind1 = (find(time > 1.782));
% ind1 = ind1(1);
% ind2 = find(time > 2.075);
% ind2 = ind2(1);
dt = time(2) - time(1);
I = trapz(time(ind1:ind2),abs(new_data(ind1:ind2)));
Isp = [Isp I/(fuel_mass*(9.81))];
delta_v = [delta_v Isp(i-2)*9.81*log(.052+fuel_mass/.052)];
Sem = zeros(1,5);
end
% hold on
% 
% plot(time,new_data)
% title('Thrust Vs Time')
% xlabel('Time(s)')
% ylabel('Thrust (N)')

% 1.652 KS/s
% end