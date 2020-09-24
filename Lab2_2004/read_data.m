function [Isp, delta_v, new_data] = read_data(filename)
data = xlsread(filename);
new_data = data(:,3);
new_data = convforce([new_data],'lbf','N');
time = linspace(0,length(data(:,3))/1652,12000);
time = time';
ind1 = (find(time > 1.782));
ind1 = ind1(1);
ind2 = find(time > 2.075);
ind2 = ind2(1);
dt = time(2) - time(1);
I = trapz(time(ind1:ind2),abs(new_data(ind1:ind2)));
Isp = I/(9.81);
delta_v = Isp*9.81*log(1.052/.052);
y = -91.5986*time + 163.1372;
hold on
plot(time,y)
plot(time,new_data)
title('Thrust Profile')
legend('Truncation Line','Thrust')
xlabel('Time (s)')
ylabel('Force (N)')

% 1.781
% 1.652 KS/s
% start @ 1.781,16.32
% done @ 2.075, -26.93
%b = 163.1372
end