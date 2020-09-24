%Author:Sage Herrin
%Created: 4/11/12
%APPM 4650 hw10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 5.5, problem 2a
clear all;clc
%begin by obtaining initial values from 4th order RK method for 2-step
%AB method
clear all;clc
t = [0:.1:1];
h = .1;
y(1) = 1;
order = 2;
%define DE
f = @(t,y) (2 - 2*t*y)/(t^2 + 1);
%run RK
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/2)*(3*f(t(j),y(j)) - f(t(j-1),y(j-1)));
end
%compare approx., exact, and error
exact = (2.*t + 1)./(t.^2 + 1)
y
error = abs(y - exact)
%repeat process for 3-stop AB method
clear y;clc
y(1) = 1;
order = 3;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/12)*(23*f(t(j),y(j)) - 16*f(t(j-1),y(j-1)) + 5*f(t(j-2),y(j-2)));
end
%compare approx., exact, and error
y
error = abs(y - exact)
%4-step order AB method
clear y;clc
y(1) = 1;
order = 4;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/24)*(55*f(t(j),y(j)) - 59*f(t(j-1),y(j-1)) + 37*f(t(j-2),y(j-2)) - 9*f(t(j-3),y(j-3)));
end
%compare approx., exact, and error
y
error = abs(y - exact)
%5-step order AB method
clear y;clc
y(1) = 1;
order = 5;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/720)*(1901*f(t(j),y(j)) - 2774*f(t(j-1),y(j-1)) + 2616*f(t(j-2),y(j-2)) - 1274*f(t(j-3),y(j-3)) + 251*f(t(j-4),y(j-4)));
end
%compare approx., exact, and error
y
error
%% Problem 2b, reusing above script and changing IVP parameters
clear all;clc
t = [1:.1:2];
h = .1;
y(1) = -log(2)^-1;
order = 2;
%define DE
f = @(t,y) y^2/(t + 1);
%run RK
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/2)*(3*f(t(j),y(j)) - f(t(j-1),y(j-1)));
end
%compare approx., exact, and error
exact = -1./(log(t + 1))
y
error = abs(y - exact)
%repeat process for 3-stop AB method
clear y;clc
y(1) = -log(2)^-1;
order = 3;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/12)*(23*f(t(j),y(j)) - 16*f(t(j-1),y(j-1)) + 5*f(t(j-2),y(j-2)));
end
%compare approx., exact, and error
y
error = abs(y - exact)
%4-step order AB method
clear y;clc
y(1) = -log(2)^-1;
order = 4;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/24)*(55*f(t(j),y(j)) - 59*f(t(j-1),y(j-1)) + 37*f(t(j-2),y(j-2)) - 9*f(t(j-3),y(j-3)));
end
%compare approx., exact, and error
y
error = abs(y - exact)
%5-step order AB method
clear y;clc
y(1) = -log(2)^-1;
order = 5;
for i = 1:(order - 1) 
    k1 = h*f(t(i),y(i));
    k2 = h*f((t(i)+h/2),(y(i)+.5*k1));
    k3 = h*f((t(i)+h/2),(y(i)+.5*k2));
    k4 = h*f(t(i+1),(y(i) + k3));
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for j = order:length(t) - 1
    y(j+1) = y(j) + (h/720)*(1901*f(t(j),y(j)) - 2774*f(t(j-1),y(j-1)) + 2616*f(t(j-2),y(j-2)) - 1274*f(t(j-3),y(j-3)) + 251*f(t(j-4),y(j-4)));
end
%compare approx., exact, and error
y
error
%% 5.9, question 2a
clear y;clc
u1(1) = -1;
u2(1) = 0;
h = 0.1;
t = [0:0.1:1];

% f1 = @(t,u1,u2) u1 - u2 + 2;
% f2 = @(t,u1,u2) -u1 +u2 +4*t;

for i = 1:length(t) - 1
    for j = 1:2
        k(1,j) = h*funcgrab1(j,t(i),u1(i),u2(i));
    end
    for j = 1:2
        k(2,j) = h*funcgrab1(j,t(i)+h/2,u1(i)+(k(1,1)*.5),u2(i)+(k(1,2)*.5));
    end
    for j = 1:2
        k(3,j) = h*funcgrab1(j,t(i)+h/2,u1(i)+(.5*k(2,1)),u2(i)+(.5*k(2,2)));
    end
    for j = 1:2
        k(4,j) = h*funcgrab1(j,t(i)+h,u1(i)+k(3,1),u2(i)+k(3,2));
    end
   u1(i+1) = u1(i) + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
   u2(i+1) = u2(i) + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;
end

u1
u2

exact_u1 = -.5*exp(2.*t) + t.^2 + 2.*t - .5
exact_u2 = .5*exp(2.*t) + t.^2 - .5
%compare approx. and actual soln. with error
error_u1 = abs(exact_u1 - u1)
error_u2 = abs(exact_u2 - u2)
%% 5.9, question 2b
clear y;clc
u1(1) = -3;
u2(1) = 5;
h = 0.2;
t = [0:0.2:2];

% f1 = @(t,u1,u2) u1 - u2 + 2;
% f2 = @(t,u1,u2) -u1 +u2 +4*t;

for i = 1:length(t) - 1
    for j = 1:2
        k(1,j) = h*funcgrab2(j,t(i),u1(i),u2(i));
    end
    for j = 1:2
        k(2,j) = h*funcgrab2(j,t(i)+h/2,u1(i)+(k(1,1)*.5),u2(i)+(k(1,2)*.5));
    end
    for j = 1:2
        k(3,j) = h*funcgrab2(j,t(i)+h/2,u1(i)+(.5*k(2,1)),u2(i)+(.5*k(2,2)));
    end
    for j = 1:2
        k(4,j) = h*funcgrab2(j,t(i)+h,u1(i)+k(3,1),u2(i)+k(3,2));
    end
   u1(i+1) = u1(i) + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
   u2(i+1) = u2(i) + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;
end

u1
u2

exact_u1 = -3*exp(t) + t.^2
exact_u2 = 4*exp(t) - 3.*t + 1
%compare approx. and actual soln. with error
error_u1 = abs(exact_u1 - u1)
error_u2 = abs(exact_u2 - u2)
%% 5.9, question 4a
clear all;clc
u1(1) = 3;
u2(1) = -1;
u3(1) = 9;
h = 0.2;
t = [0:0.2:2];

% f1 = @(t,u1,u2) u1 - u2 + 2;
% f2 = @(t,u1,u2) -u1 +u2 +4*t;

for i = 1:length(t) - 1
    for j = 1:3
        k(1,j) = h*funcgrab3(j,t(i),u1(i),u2(i),u3(i));
    end
    for j = 1:3
        k(2,j) = h*funcgrab3(j,t(i)+h/2,u1(i)+(k(1,1)*.5),u2(i)+(k(1,2)*.5),u3(i)+(k(1,3)*.5));
    end
    for j = 1:3
        k(3,j) = h*funcgrab3(j,t(i)+h/2,u1(i)+(.5*k(2,1)),u2(i)+(.5*k(2,2)),u3(i)+(.5*k(2,3)));
    end
    for j = 1:3
        k(4,j) = h*funcgrab3(j,t(i)+h,u1(i)+k(3,1),u2(i)+k(3,2),u3(i)+k(3,3));
    end
   u1(i+1) = u1(i) + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
   u2(i+1) = u2(i) + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;
   u3(i+1) = u3(i) + (k(1,3) + 2*k(2,3) + 2*k(3,3) + k(4,3))/6;
end
u1
exact = exp(-t) + exp(2.*t) + exp(-2.*t)
error = abs(exact - u1)
%% 5.9, question 4b
clear all;clc
u1(1) = 2;
u2(1) = 8;
u3(1) = 6;
h = 0.1;
t = [1:0.1:2];

% f1 = @(t,u1,u2) u1 - u2 + 2;
% f2 = @(t,u1,u2) -u1 +u2 +4*t;

for i = 1:length(t) - 1
    for j = 1:3
        k(1,j) = h*funcgrab4(j,t(i),u1(i),u2(i),u3(i));
    end
    for j = 1:3
        k(2,j) = h*funcgrab4(j,t(i)+h/2,u1(i)+(k(1,1)*.5),u2(i)+(k(1,2)*.5),u3(i)+(k(1,3)*.5));
    end
    for j = 1:3
        k(3,j) = h*funcgrab4(j,t(i)+h/2,u1(i)+(.5*k(2,1)),u2(i)+(.5*k(2,2)),u3(i)+(.5*k(2,3)));
    end
    for j = 1:3
        k(4,j) = h*funcgrab4(j,t(i)+h,u1(i)+k(3,1),u2(i)+k(3,2),u3(i)+k(3,3));
    end
   u1(i+1) = u1(i) + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6;
   u2(i+1) = u2(i) + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6;
   u3(i+1) = u3(i) + (k(1,3) + 2*k(2,3) + 2*k(3,3) + k(4,3))/6;
end
u1
exact = 2.*t - t.^-1 + t.^2 + t.^3 - 1
error = abs(exact - u1)
        
        






