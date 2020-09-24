clear all;close all;clc

C_l = 1.3518;
mars_rho = 0.02; %Kg/m^3
mars_weight = 6.4*3.71; %Newtons
C_D0 = 0.0246;
e = 0.6010;
AR = 16.5;

for i=1:20
S=.63*[.5:.5:10];
S=S(i);
P_a = 26.1*50;
mars_weight = 6.4*3.71;
syms v_max
v_stall(i)=sqrt(2*mars_weight/(max(C_l)*mars_rho*S));
P_r =.5*mars_rho*(v_max)^3*S*C_D0 + ((mars_weight)^2/(.5*mars_rho*v_max*S*pi*e*AR));
b=sqrt(AR*S);
store=double(solve(P_r == P_a,v_max));
A(i) = (store(2))
end