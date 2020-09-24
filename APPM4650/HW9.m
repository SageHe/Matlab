%problem 2d of 5.3
clear all;clc
t = [1:.25:1.75];
y(1) = 2;
h = .25;
for i = 1:length(t)
    y(i+1) = y(i) + .25*(t(i)^-2*(sin(2*t(i)) - 2*t(i)*y(i)) + .125*(-2*t(i)^-3*(sin(2*t(i)) - 2*t(i)*y(i))...
        + t(i)^-2*(2*cos(2*t(i)) - 2*y(i) - 2*t(i)*(t(i)^-2*sin(2*t(i)) - 2*t(i)*y(i)))));
end
%problem 9a of 5.3
clear all;clc
t = [1:.1:1.9];
y1(1) = 0;
h = .1;
for i = 1:length(t)
    y1(i+1) = y1(i) + .1*((2/t(i))*y1(i) + t(i)^2*exp(t(i)) + .05*((-2/t(i)^2)*y1(i) + 2/t(i)*((2/t(i))*y1(i) + t(i)^2*exp(t(i))) + 2*t(i)*exp(t(i)) + t(i)^2*exp(t(i))));
end
%evaluation of exact function
t = [1:.1:2];
y2 = t.^2.*(exp(t) - exp(1));
%problem 14a of 5.4
clear all;clc
t = [0:.5:1];
h = .5;
y(1) = 1;
for i = 1:length(t) - 1
    k1 = h*exp(t(i) - y(i));
    k2 = h*exp((t(i)+h/2) - (y(i)+.5*k1));
    k3 = h*exp((t(i)+h/2) - (y(i)+.5*k2));
    k4 = h*exp(t(i+1) - (y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
%problem 14b of 5.4
clear all;clc
t = [1:.5:2];
h = .5;
y(1) = 2;
for i = 1:length(t) - 1
    k1 = h*((1+t(i))/(1+y(i)));
    k2 = h*((1+t(i)+.5*h)/(1+y(i)+.5*k1));
    k3 = h*((1+t(i)+.5*h)/(1+y(i)+.5*k2));
    k4 = h*((1+t(i+1))/(1+y(i)+k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
    
    
    
    





