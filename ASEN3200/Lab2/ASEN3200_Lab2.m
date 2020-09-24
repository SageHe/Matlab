data = [40 7.12;30 6.08;25 4.65;20 3.81]; %first column velocity, second time period
spin_rate = (1./data(:,2))*2*pi; %Determine precession rate from period
w_x = data(:,1)*1000*(1/3600); %determine angular vel. from measured linear vel.
w_x = w_x./.3302;
% figure
% grid on
% grid minor
    