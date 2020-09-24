clear all;close all;clc
r2 = [30 40 50 75 100 112.5 125];
r2 = r2/1000;
r1 = r2 - .01;
L = [1 2 3 4 5 6 7];
I = pi/4*(r2.^4-r1.^4);
E = 200e9;
P_cr = (pi^2.*E.*I)./L.^2;

A = pi.*(r2.^2 - r1.^2);
sig = P_cr./A;

sig_y = 300e6;
short_long = [];
for i = 1:numel(sig)
    if sig(i)>sig_y
        fprintf('Bar %d is short\n',i)
    else
        fprintf('Bar %d is long\n',i)
    end
end

        