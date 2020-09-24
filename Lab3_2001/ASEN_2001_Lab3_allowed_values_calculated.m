clear; clc; format compact; close all;

data = xlsread('TestData_edited.xlsx');
num = 18;
test_no = data(:,1);
F = data(:,2); % Force
a = data(:,3); % distance from left most support 
w = data(:,4); % width of beam top down view
d_f = data(:,5); % distance from the failure location of beam to closet support
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_b = 2*(((w.*.00079375^3)/12)+(w.*.00079375*.009906^2));
I_f = (w.*0.01905^3)/12;

E_b = 3.2953*10^9;
E_f = 0.035483 * 10^9;
I = (I_b + (E_f/E_b)*I_f)*0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0.0216375; % Thickness [m]
c = T/2; % Half the cross section's total height

Area = w*T; % Area [m^2]

tau_n = 3/4*F./Area; % Shear stress [pa]

M_b_a = (.5*F).*a; % Bending moment
M_b_d_f = (.5*F).*d_f;

for i=1:num
   if a(i) < d_f(i)
       M_b(i,1) = M_b_a(i);
   else
       M_b(i,1) = M_b_d_f(i);
   end
end

sig_n = (M_b.*c)./I; % Bending stress
break_shear = zeros(1,1);
j=1;
k=1;
INDEX = zeros(num,1);
for i=1:num
    if d_f(i) <= a(i)
            %broke in shear
           INDEX(i,1) = 9999999; 
        break_shear(j) = tau_n(i);
        j=j+1;
    else 
        break_mom(k) = sig_n(i);
        INDEX(i,1) = 6666666666;
        k=k+1;
    end
end
X_tau = [1:7];
X_sig = [1:11];
tau_fit = polyfit(X_tau,break_shear,1);
sig_fit = polyfit(X_sig,break_mom,1);
plot(break_shear,'*')
hold on
tau_fit_eqn = X_tau * tau_fit(1)+tau_fit(2);
plot(X_tau,tau_fit_eqn);
title('shear')
figure
plot(break_mom,'*')
hold on
sigma_fit = X_sig * sig_fit(1)+sig_fit(2);
plot(X_sig,sigma_fit)
title('moment')

figure 
plot(I,'*')

tau_avg = mean(break_shear);
sig_avg = mean(break_mom);
I_avg = mean(I);
I_b_avg = mean(I_b);
I_f_avg = mean(I_f);
tau_allow = tau_avg/1.5;
sig_allow = sig_avg/1.5;


distance = 0.9144;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I_b = 3.35*10^-11;
% I_f = 5.85*10^-8;
% I = I_b + (E_f/E_b)*I_f;
% I_avg = mean(I);
%sig_allow = 3.58*10^8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sig_allow = 3.58*10^8;
syms p0 x L m
p = p0 * sqrt(1 - (2 * x / distance)^2);
reacts = -int(p,x,-0.5*distance,0.5*distance) * 0.5;
shearF = int(p,x,-0.5*distance,x) + reacts;
momentF = int(shearF,x,-0.5*distance,x);
momentF = subs(momentF,x,0);
momentF_mult = double(subs(momentF,p0,1));
P0 = -(sig_allow*I_avg)/(c*momentF_mult*.1016);

