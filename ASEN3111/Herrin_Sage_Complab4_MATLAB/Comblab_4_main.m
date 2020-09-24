%Problem 1 , calculate e,c_L, c_Di from inputs to function
clear all;close all; clc
set(0,'defaulttextInterpreter','latex');
[e, c_L, c_Di] = PLLT(100,7.0459,7.0326,5,15,0,-.0370,0,.0873,1000);
%Problem 2, calculate the lift and induced drag on the wing at a sea level
%velocity of 150 mph and calculate number of odd terms needed for
%5%,1%, and .1% error
c_r = 15;
c_t = 5;
b = 100;
S = ((c_r + c_t)/2)*b;
v_inf = 150;
v_inf = convvel(150,'mph','ft/s');
rho = .002378; %slug/ft^3
q_inf = .5*rho*v_inf^2;
L = c_L*S*q_inf;
D_i = c_Di*S*q_inf;
fprintf('The lift for the wing at a velocity of 150 miles/hour at sea level is %.3f pounds\n',L)
fprintf('The induced drag for the wing at a velocity of 150 miles/hour at sea level is %.3f pounds\n',D_i)
accepted_L = L;
accepted_D = D_i;
%calculate number of odd terms for 5% relative error
[e, c_L, c_Di] = PLLT(100,7.0459,7.0326,5,15,0,-.0370,0,.0873,1);
L = q_inf*S*c_L;
D_i = q_inf*S*c_Di;
i = 2;
L_err_05 = [];
L_err_01 = [];
D_err_05 = [];
D_err_01 = [];
while abs((accepted_L - L(end))/accepted_L) > .001
    [e, c_L, c_Di] = PLLT(100,7.0459,7.0326,5,15,0,-.0370,0,.0873,i);
    L = [L q_inf*S*c_L];
    if abs((accepted_L - L(end)))/accepted_L < .05 && isempty(L_err_05)
        L_err_05 = i;
    end
    if abs((accepted_L - L(end)))/accepted_L < .01 && isempty(L_err_01)
        L_err_01 = i;
    end
    i = i + 1;
end
L_err_001 = i;
i = 2;
while abs((accepted_D - D_i(end))/accepted_D) > .001
    [e, c_L, c_Di] = PLLT(100,7.0459,7.0326,5,15,0,-.0370,0,.0873,i);
    D_i = [D_i q_inf*S*c_Di];
    if abs((accepted_D - D_i(end)))/accepted_D < .05 && isempty(D_err_05)
        D_err_05 = i;
    end
    if abs((accepted_D - D_i(end)))/accepted_D < .01 && isempty(D_err_01)
        D_err_01 = i;
    end
    i = i + 1;
end
D_err_001 = i;
fprintf('It takes %d odd terms for lift to converge within 5%% relative error\n',L_err_05)
fprintf('It takes %d odd terms for lift to converge within 1%% relative error\n',L_err_01)
fprintf('It takes %d odd terms for lift to converge within .1%% relative error\n',L_err_001)
fprintf('It takes %d odd terms for drag to converge within 5%% relative error\n',D_err_05)
fprintf('It takes %d odd terms for drag to converge within 1%% relative error\n',D_err_01)
fprintf('It takes %d odd terms for drag to converge within .1%% relative error\n',D_err_001)

figure(1)
hold on
grid on
% grid minor
grid minor
plot([1:numel(L)],L)
plot([1:numel(L)],ones(1,numel(L))*1.05*accepted_L)
plot([1:numel(L)],ones(1,numel(L))*1.01*accepted_L)
plot([1:numel(L)],ones(1,numel(L))*1.001*accepted_L)
set(gca,'TickLabelInterpreter','latex');
title('Relative Error in Lift VS Number of Odd Terms')
xlabel('Number of Odd Terms')
ylabel('Lift [pounds]')
legend('Running Lift Value','5%','1%','.1%')
lgd = legend;
lgd.FontSize = 7;
lgd.Title.String = 'Relative Error in Lift';

figure(2)
hold on
grid on
grid minor
plot([1:numel(D_i)],D_i)
plot([1:numel(D_i)],ones(1,numel(D_i))*1.05*accepted_D)
plot([1:numel(D_i)],ones(1,numel(D_i))*1.01*accepted_D)
plot([1:numel(D_i)],ones(1,numel(D_i))*1.001*accepted_D)
set(gca,'TickLabelInterpreter','latex');
title('Relative Error in Drag VS Number of Odd Terms')
xlabel('Number of Odd Terms')
ylabel('Drag [pounds]')
legend('Running Drag Value','5%','1%','.1%')
lgd = legend;
lgd.FontSize = 7;
lgd.Title.String = 'Relative Error in Drag';


%Calculate and plot e vs ct/cr for given AR's at fixed geometric angle of
%attack of 3 degrees
AR = [4 6 8 10];
c_r = 15;
c_t = linspace(.001,15,1000);
for i = 1:length(AR)
    for j = 1:length(c_t)
        b = AR(i)*((c_t(j) + c_r)/2);
        [e, c_L, c_Di] = PLLT(b,2*pi,2*pi,c_t(j),15,0,0,0.0524,0.0524,50);
        e_vec(i,j) = e;
    end
end
figure(3)
hold on
grid on
grid minor 
plot((c_t./c_r),e_vec(1,:))
plot((c_t./c_r),e_vec(2,:))
plot((c_t./c_r),e_vec(3,:))
plot((c_t./c_r),e_vec(4,:))
set(gca,'TickLabelInterpreter','latex');
title('Span Efficfiency Factor VS Taper Ratio')
xlabel('Taper Ratio')
ylabel('Span Efficiency Factor e')
legend('4','6','8','10','location','southeast')
lgd = legend;
lgd.FontSize = 7;
lgd.Title.String = 'Aspect Ratio';

delta = (1 - e_vec)./e_vec;

figure(4)
hold on
grid on
grid minor
plot((c_t./c_r),delta(1,:))
plot((c_t./c_r),delta(2,:))
plot((c_t./c_r),delta(3,:))
plot((c_t./c_r),delta(4,:))
set(gca,'TickLabelInterpreter','latex');
title('delta VS taper ratio')
xlabel('Taper Ratio')
ylabel('delta')
legend('4','6','8','10','location','northeast')
lgd = legend;
lgd.FontSize = 7;
lgd.Title.String = 'Aspect Ratio';
