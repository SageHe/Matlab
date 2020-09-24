clear all;
close all;
clc;

[t,y]= ode45('react',[0 4],[1 0 0]);

plot(t,y)
