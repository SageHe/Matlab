clear; clc; format compact; close all;

data = xlsread('TestData.xlsx');

test_no = data(:,1);
F = data(:,2); % Force
a = data(:,3); % distance from left most support 
w = data(:,4); % width of beam top down view
d_f = data(:,5); % distance from the failure location of beam to closet support
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_b = 2*(((w.*.00079375^3)/12)+w.*.00079375*.00992^2);
I_f = (w.*0.01905^3)/12;

E_b = 3.2953*10^9;
E_f = 0.035483 * 10^9;
I = I_b + (E_f/E_b)*I_f;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0.0206375; % Thickness [m]
c = T/2; % Half the cross section's total height

Area = w*T; % Area [m^2]

tau_n = (3/4)*F./Area; % Shear stress [pa]

M_b_a = (.5*F).*a; % Bending moment
M_b_d_f = (.5*F).*d_f;

for i=1:length(Area)
   if a(i) < d_f(i)
       M_b(i,1) = M_b_a(i);
   else
       M_b(i,1) = M_b_d_f(i);
   end
end

sig_n = (M_b.*c)./I; % Bending stress

for i=1:length(sig_n)
    if d_f(i) < a(i)
        if tau_n(i) < sig_n(i)
            %broke in moment
            break_mom(i) = sig_n(i);
        else
            %broke in shear
            break_shear(i) = tau_n(i); 
        end
    else 
        break_mom(i) = sig_n(i);
    end
end

