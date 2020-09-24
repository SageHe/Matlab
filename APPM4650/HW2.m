clear all;close all;clc
gamma = 0.5772156649;
gamma_n = 0;
g_n = 0;
n = 1;
while floor(g_n*10^4)/10^4 ~= .5772
    gamma_n = [gamma_n (1/n)];
    g_n = sum(gamma_n) - log(n);
    n = n + 1;
end

n = n - 1;