%Load original data from 'data' file, parse out incorrect data files and
%narrow down to single repition for airfoil and cylinder, save to .mat
%file, comment out original read-in process and parse data for speed
tic
clear all;close all;clc
% Full_Data = dir('Data');
% 
% %Find third repetition of downstream airfoil data 
% j = 0;
% for i = 3:size(Full_Data,1)
%     if isempty(strfind(Full_Data(i).name,'r03')) == 0 && isempty(strfind(Full_Data(i).name,'Airfoil')) == 0 && isempty(strfind(Full_Data(i).name,'Down')) == 0
%         j = j + 1;
%         airfoil_data_down(j).name = Full_Data(i).name;
%         temp_data = load(airfoil_data_down(j).name);
%         temp_data = temp_data(:,[3 4 6 27 28]);
%         airfoil_data_down(j).data = temp_data;
%     end
% end
% %Find third repetition of upstream airfoil data 
% j = 0;
% for i = 3:size(Full_Data,1)
%     if isempty(strfind(Full_Data(i).name,'r03')) == 0 && isempty(strfind(Full_Data(i).name,'Airfoil')) == 0 && isempty(strfind(Full_Data(i).name,'Up')) == 0
%         j = j + 1;
%         airfoil_data_up(j).name = Full_Data(i).name;
%         temp_data = load(airfoil_data_up(j).name);
%         temp_data = temp_data(:,[3 4 6 27 28]);
%         airfoil_data_up(j).data = temp_data;
%     end
% end
% %Find fourth repetition of upstream cylinder data
% j = 0;
% for i = 3:size(Full_Data,1)
%     if isempty(strfind(Full_Data(i).name,'r04')) == 0 && isempty(strfind(Full_Data(i).name,'Cylinder')) == 0 && isempty(strfind(Full_Data(i).name,'Up')) == 0
%         j = j + 1;
%         cylinder_data_up(j).name = Full_Data(i).name;
%         temp_data = load(cylinder_data_up(j).name);
%         temp_data = temp_data(:,[3 4 6 27 28]);
%         cylinder_data_up(j).data = temp_data;
%     end
% end
% %Find fourth repetition of downstream cylinder data
% j = 0;
% for i = 3:size(Full_Data,1)
%     if isempty(strfind(Full_Data(i).name,'r04')) == 0 && isempty(strfind(Full_Data(i).name,'Cylinder')) == 0 && isempty(strfind(Full_Data(i).name,'Down')) == 0
%         j = j + 1;
%         cylinder_data_down(j).name = Full_Data(i).name;
%         temp_data = load(cylinder_data_down(j).name);
%         temp_data = temp_data(:,[3 4 6 27 28]);
%         cylinder_data_down(j).data = temp_data;
%     end
% end
% % cylinder_data_down(6).name = Full_Data(125).name;
% % cylinder_data_down(6).data = load(Full_Data(125).name);
% cylinder_data_down(6).name = Full_Data(124).name;
% cylinder_data_down(6).data = load(Full_Data(124).name);
% cylinder_data_down(13).name = Full_Data(152).name;
% cylinder_data_down(13).data = load(Full_Data(152).name);
% cylinder_data_down(6).data = cylinder_data_down(6).data(:,[3 4 6 27 28]);
% cylinder_data_down(13).data = cylinder_data_down(13).data(:,[3 4 6 27 28]);
% for i = 1:size(airfoil_data_up,2)
%     temp_data_airfoil_up = airfoil_data_up(i).data;
%     temp_data_airfoil_down = airfoil_data_down(i).data;
%     temp_data_cylinder_up = cylinder_data_up(i).data;
%     temp_data_cylinder_down = cylinder_data_down(i).data;
%     pos = [1 500 501 1000 1001 1500 1501 2000 2001 2500 2501 3000 3001 3500 3501 4000 4001 4500 4501 5000 5001 5500 5501 6000 6001 6500 6501 7000 7001 7500 7501 8000 8001 8500 8501 9000 9001 9500 9501 10000];
%     new_temp_airfoil_data_up = zeros(20,5);
%     new_temp_airfoil_data_down = zeros(20,5);
%     new_temp_cylinder_data_up = zeros(20,5);
%     new_temp_cylinder_data_down = zeros(20,5);
%     for j = 1:(length(pos)/2)
%         new_temp_airfoil_data_up(j,:) = mean(temp_data_airfoil_up([pos(2*j-1):pos(2*j)],:),1);
%         new_temp_airfoil_data_down(j,:) = mean(temp_data_airfoil_down([pos(2*j-1):pos(2*j)],:),1);
%         new_temp_cylinder_data_up(j,:) = mean(temp_data_cylinder_up([pos(2*j-1):pos(2*j)],:),1);
%         new_temp_cylinder_data_down(j,:) = mean(temp_data_cylinder_down([pos(2*j-1):pos(2*j)],:),1);
%     end
%     airfoil_data_up(i).data = new_temp_airfoil_data_up;
%     airfoil_data_down(i).data = new_temp_airfoil_data_down;
%     cylinder_data_up(i).data = new_temp_cylinder_data_up;
%     cylinder_data_down(i).data = new_temp_cylinder_data_down;
% end
% airfoil_data_down(6).data(:,4) = 38;
load('pert_data');
%Downstream airfoil @ 15 m/s
max_vel_def = [];
x = [];
hww_airfoil_15 = [];
Cd = [];
C = 8.89*10^(-2); %Airfoil chord length in meters
for h = 1:(size(airfoil_data_down,2)/2)
    data = airfoil_data_down(h).data;
    vel_def = [];
    for i = 1:size(data,1)
        q = data(i,3);
        rho = data(i,1);
        v = sqrt((2*q)/rho);
        vel_def = [vel_def (data(i,2) - v)];
    end
    y = data(:,5);
    q_avg = data(:,3);
    rho_avg = data(:,1);
    v_avg = sqrt((2.*q_avg)./rho_avg);
    V_avg = mean(v_avg); %Use this for actual drag calc
    D = trapz(y(:)*10^(-3),(V_avg*vel_def*(mean(data(:,1)))));
    Cd = [Cd (D/((.5*mean(data(:,1))*(V_avg)^2)*C))];
    figure(1)
    hold on
   subplot(2,2,1), plot(data(:,5),vel_def)
    hold off
    title('Velocity Deficit of Airfoil @ 15 m/s')
    xlabel('Vertical Location (mm)')
    ylabel('Velocity Deficit (m/s)')
    legend('x = 13mm','x = 18mm','x = 23mm','x = 28mm','x = 33mm','x = 38mm',...
    'x = 43mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Airfoil V = 15';
    %Calc max velocity deficit
    max_vel_def = [max_vel_def max(vel_def)];
    x = [x mean(data(:,4));];
    %Calc half wake width
    above = find(vel_def > (max_vel_def(end))/2);
    A = above(end) + 1;
    B = above(1) - 1;
    hww_airfoil_15 = [hww_airfoil_15 .5*(data(A,5) - data(B,5))];
    %Calc non-dimensional velocity deficit
    nd_vel_def = (vel_def./max_vel_def(end));
    nd_ver_pos = y/hww_airfoil_15(end);
    figure(3)
    hold on
    subplot(2,2,1), plot(nd_ver_pos,nd_vel_def)
    title('Non-dimensional Velocity VS Non-dimensional Vertical Location, Airfoil @ 15 m/s')
    xlabel('Non-dimensional Vertical Location')
    ylabel('Non-dimensional Velocity')
    legend('x = 13mm','x = 18mm','x = 23mm','x = 28mm','x = 33mm','x = 38mm',...
    'x = 43mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Airfoil V = 15';
end
Cd_Airfoil_15 = mean(Cd);
figure(2)
subplot(2,2,1), plot(x,max_vel_def,'o')
title('Max velocity Deficit for Airfoil @ 15 m/2')
xlabel('X Position (mm')
ylabel('Velocity Deficit (m/s)')
figure(4)
hold on
subplot(2,2,1), plot(x,hww_airfoil_15)
title('Half Wake Width VS Horizontal Position, Airfoil @ 15 m/s')
xlabel('Horizontal Position (mm)')
ylabel('Half Wake Width (mm)')
%Downstream airfoil @ 25 m/s
max_vel_def = [];
x = [];
hww_airfoil_25 = [];
Cd = [];
for h = ((size(airfoil_data_down,2)/2) + 1):size(airfoil_data_down,2)
    data = airfoil_data_down(h).data;
    vel_def = [];
    for i = 1:size(data,1)
        q = data(i,3);
        rho = data(i,1);
        v = sqrt((2*q)/rho);
        vel_def = [vel_def (data(i,2) - v)];
    end
    y = data(:,5);
    q_avg = data(:,3);
    rho_avg = data(:,1);
    v_avg = sqrt((2.*q_avg)./rho_avg);
    V_avg = mean(v_avg); %Use this for actual drag calc
    D = trapz(y(:)*10^(-3),(V_avg*vel_def*(mean(data(:,1)))));
    Cd = [Cd (D/((.5*mean(data(:,1))*(V_avg)^2)*C))];
    figure(1)
    hold on
    subplot(2,2,2), plot(data(:,5),vel_def)
    hold off
        title('Velocity Deficit of Airfoil @ 25 m/s')
    xlabel('Vertical Location (mm)')
    ylabel('Velocity Deficit (m/s)')
    legend('x = 13mm','x = 18mm','x = 23mm','x = 28mm','x = 33mm','x = 38mm',...
    'x = 43mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Airfoil V = 25';
    max_vel_def = [max_vel_def max(vel_def)];
    x = [x mean(data(:,4));];
    y = data(:,5);
    above = find(vel_def > (max_vel_def(end))/2);
    A = above(end) + 1;
    B = above(1) - 1;
    hww_airfoil_25 = [hww_airfoil_25 .5*(data(A,5) - data(B,5))];
    nd_vel_def = (vel_def./max_vel_def(end));
    nd_ver_pos = y/hww_airfoil_25(end);
    figure(3)
    hold on
    subplot(2,2,2), plot(nd_ver_pos,nd_vel_def)
    title('Non-dimensional Velocity VS Non-dimensional Vertical Position, Airfoil @ 25 m/s')
    xlabel('Non-dimensional Vertical Position')
    ylabel('Non-dimensional Velocity')
    legend('x = 13mm','x = 18mm','x = 23mm','x = 28mm','x = 33mm','x = 38mm',...
    'x = 43mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Airfoil V = 25';

end
Cd_airfoil_25 = mean(Cd);
figure(2)
subplot(2,2,2), plot(x,max_vel_def,'o')
title('Max Velocity Deficit of Airfoil @ 25 m/s')
xlabel('X Position (mm)')
ylabel('Velocity Deficit')
figure(4)
hold on
subplot(2,2,2), plot(x,hww_airfoil_25)
title('Half Wake Width VS Horizontal Position, Airfoil @ 25 m/s')
xlabel('Horizontal Position (mm)')
ylabel('Half Wake Width (mm)')
%Downstream cylinder @ 15 m/s
max_vel_def = [];
x = [];
hww_cylinder_15 = [];
Cd = [];
C = 1.27*10^(-2);
for h = 1:(size(cylinder_data_down,2)/2)
    data = cylinder_data_down(h).data;
    vel_def = [];
    for i = 1:size(data,1)
        q = data(i,3);
        rho = data(i,1);
        v = sqrt((2*q)/rho);
        vel_def = [vel_def (data(i,2) - v)];
    end
    y = data(:,5);
    q_avg = data(:,3);
    rho_avg = data(:,1);
    v_avg = sqrt((2.*q_avg)./rho_avg);
    V_avg = mean(v_avg); %Use this for actual drag calc
    D = trapz(y(:)*10^(-3),(V_avg*vel_def*(mean(data(:,1)))));
    Cd = [Cd (D/((.5*mean(data(:,1))*(V_avg)^2)*C))];
    figure(1)
    hold on
    subplot(2,2,3), plot(data(:,5),vel_def)
    hold off
    title('Velocity Deficit of cylinder @ 15 m/s')
    xlabel('Vertical Location (mm)')
    ylabel('Velocity Deficit (m/s)')
    legend('x = 120mm','x = 150mm','x = 180mm','x = 210mm','x = 240mm','x = 60mm',...
    'x = 90mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Cylinder V = 15';
    max_vel_def = [max_vel_def max(vel_def)];
    x = [x mean(data(:,4));];
    y = data(:,5);
    above = find(vel_def > (max_vel_def(end))/2);
    A = above(end) + 1;
    B = above(1) - 1;
    hww_cylinder_15 = [hww_cylinder_15 .5*(data(A,5) - data(B,5))];
    nd_vel_def = (vel_def./max_vel_def(end));
    nd_ver_pos = y/hww_cylinder_15(end);
    figure(3)
    hold on
    subplot(2,2,3), plot(nd_ver_pos,nd_vel_def)
    title('Non-dimensional Velocity VS Non-dimensional Vertical Position, Cylinder @ 15 m/s')
    xlabel('Non-dimensional Vertical Position')
    ylabel('Non-Dimensional Velocity')
    legend('x = 60mm','x = 90mm','x = 120mm','x = 150mm','x = 180mm','x = 210mm',...
    'x = 240mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Cylinder V = 15';
end
Cd_cylinder_15 = mean(Cd);
x = x([6 7 1 2 3 4 5]);
max_vel_def = max_vel_def([ 6 7 1 2 3 4 5]);
hww_cylinder_15 = hww_cylinder_15([6 7 1 2 3 4 5]);
figure(2)
subplot(2,2,3), plot(x,max_vel_def,'o')
title('Velocity Deficit of Cylinder @ 15 m/s')
xlabel('X Position (mm)')
ylabel('Velocity Deficit (m/s)')
figure(4)
hold on
subplot(2,2,3), plot(x,hww_cylinder_15)
title('Half Wake Width VS Horizontal Position, Cylinder @ 15 m/s')
xlabel('Horizontal Position (mm)')
ylabel('Half Wake Width (mm)')
%Donwstream cylinder @ 25 m/s
max_vel_def = [];
x = [];
hww_cylinder_25 = [];
Cd = [];
for h = ((size(cylinder_data_down,2)/2) + 1):size(cylinder_data_down,2)
    data = cylinder_data_down(h).data;
    vel_def = [];
    for i = 1:size(data,1)
        q = data(i,3);
        rho = data(i,1);
        v = sqrt((2*q)/rho);
        vel_def = [vel_def (data(i,2) - v)];
    end
    y = data(:,5);
    q_avg = data(:,3);
    rho_avg = data(:,1);
    v_avg = sqrt((2.*q_avg)./rho_avg);
    V_avg = mean(v_avg); %Use this for actual drag calc
    D = trapz(y(:)*10^(-3),(V_avg*vel_def*(mean(data(:,1)))));
    Cd = [Cd (D/((.5*mean(data(:,1))*(V_avg)^2)*C))];
    figure(1)
    hold on
    subplot(2,2,4), plot(data(:,5),vel_def)
    hold off
    title('Velocity Deficit of Cylinder @ 25 m/s')
    xlabel('Vertical Location (mm)')
    ylabel('Velocity Deficit (m/s)')
    legend('x = 120mm','x = 150mm','x = 180mm','x = 210mm','x = 240mm','x = 60mm',...
    'x = 90mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Cylinder V = 25';
    max_vel_def = [max_vel_def max(vel_def)];
    x = [x mean(data(:,4));];
    y = data(:,5);
    above = find(vel_def > (max_vel_def(end))/2);
    A = above(end) + 1;
    B = above(1) - 1;
    hww_cylinder_25 = [hww_cylinder_25 .5*(data(A,5) - data(B,5))];
    nd_vel_def = (vel_def./max_vel_def(end));
    nd_ver_pos = y/hww_cylinder_25(end);
    figure(3)
    hold on
    subplot(2,2,4), plot(nd_ver_pos,nd_vel_def)
    title('Non-dimensional Velocity VS Non-dimensional Vertical Position, Cylinder @ 25 m/s')
    xlabel('Non-dimensional Vertical Position')
    ylabel('No-dimensional Velocity')
    legend('x = 60mm','x = 90mm','x = 120mm','x = 150mm','x = 180mm','x = 210mm',...
    'x = 240mm')
    lgd = legend;
    lgd.FontSize = 7;
    lgd.Title.String = 'Airfoil V = 15';

end
Cd_cylinder_25 = mean(Cd);
x = x([6 7 1 2 3 4 5]);
max_vel_def = max_vel_def([6 7 1 2 3 4 5]);
hww_cylinder_25 = hww_cylinder_25([6 7 1 2 3 4 5]);
figure(2)
subplot(2,2,4), plot(x,max_vel_def,'o')
title('Velocity Deficit of Cylinder @ 25 m/s')
xlabel('X Position (mm)')
ylabel('Velocity Deficit (m/s)')
figure(4)
hold on
subplot(2,2,4), plot(x,hww_cylinder_25)
title('Half Wake Width VS Horizontal Position, Cylinder @ 25 m/s')
xlabel('Horizontal Position (mm)')
ylabel('Half Wake Width (mm)')
toc
    