function [out] = calc_gamma(x,c,alpha,V_inf)
num = 1 - (x/c);
denom = x/c;
out = 2*alpha*V_inf*sqrt(num/denom);