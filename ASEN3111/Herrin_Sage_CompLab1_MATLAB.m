%Sage Herrin
%106071909
%ASEN 3111, Aerodynamics
%Comp Lab 1
%Created on 9/7/2018

%Housekeeping
clear all;close all;clc
%%Part one
%Define given constants of freesteam pressure, density, and velocity
rho = 1.225; %Kg/m^3
p_inf = 101.3e3; %Pa
V_inf = 30; %m/s
%Calculate dynamics pressure
q = .5*rho*(V_inf)^2; %Dynamic pressure
%Given two equation for C_p, p can be solved for in terms of theta
N = 100; %# of subintervals
theta = linspace(0,2*pi,(N + 1));
p = (q*(1 - 4.*(sin(theta)).^2)) + p_inf;
%Calculate N' and A' using previously calculated values of p across theta
%range 
R = .5; %Radius of circle
N_prime = p.*sin(theta)*R;
A_prime = p.*cos(theta)*R;
%Integrate N_prime and A_prime to determine lift and drag using Simpson's
%rule
sphere_lift = 0; %Preallocate lift
sphere_drag = 0;%Preallocate drag
for i = 1:(N/2)
    sphere_lift = (sphere_lift + N_prime(2*i - 1) + 4*N_prime(2*i) + N_prime(2*i + 1));
    sphere_drag = (sphere_drag + A_prime(2*i - 1) + 4*A_prime(2*i) + A_prime(2*i + 1));
end
h = (2*pi)/N;
sphere_lift = sphere_lift*((h/3)*R);
sphere_drag = sphere_drag*((h/3)*R);

%%Question 2
%Load given spline variables
load Cp
%Define airfoil parameters
C = 2; %meters
alpha = 9; %degrees
%# of subintervals
N = 50000;
x = linspace(0,C,N + 1);
yt = (.12/.2) .* C .* (.2969.*sqrt(x./C) - .1260.*(x./C) - ...
.3516.*(x./C).^2 + .2843.*(x./C).^3 - .1036.*(x./C).^4);
% hold on
% plot(x,yt)
% plot(x,-yt)
% axis equal
%Determine upper and lower pressure values
cp_upper = fnval(Cp_upper,x/C);
cp_lower = fnval(Cp_lower,x/C);
p_upper = cp_upper*q + p_inf;
p_lower = cp_lower*q + p_inf;
N_prime = 0;
A_prime = 0;
for j = 1:N 
    N_prime = (N_prime + (-((p_upper(j + 1) + p_upper(j))/2)*((x(j + 1) - x(j)))) + ((p_lower(j + 1) + p_lower(j))/2)*(x(j + 1) - x(j)));
    A_prime = (A_prime + (((p_upper(j + 1) + p_upper(j))/2)*(yt(j + 1) - yt(j))) + ((p_lower(j + 1) + p_lower(j))/2)*(yt(j + 1) - yt(j)));
end
Lift = N_prime*cosd(alpha) - A_prime*sind(alpha);
Drag = N_prime*sind(alpha) + A_prime*cosd(alpha);

fprintf('Lift is calculated to %.4f Newtons of drag\n',Lift);
fprintf('Drag is calculated to %.4f Newtons of drag\n',Drag);

%Running the above script with an N value of 10,000,000 yields results of
%lift = 1188.1095 N and drag = 1.4655 N which will be used to compare for
%relative error

N_prime = 0;
A_prime = 0;

Lift = N_prime*cosd(alpha) - A_prime*sind(alpha);
Drag = N_prime*sind(alpha) + A_prime*cosd(alpha);

N = 20;

accepted_lift = 1188.1095; %Newtons
accepted_drag = 1.4655; %Netwons

while ((Lift <= 0) || (((accepted_lift - Lift)/accepted_lift) > .05) && ((accepted_lift - Lift)/accepted_lift) > 0)
    N = N + 1;
    x = linspace(0,C,N + 1);
    yt = (.12/.2) .* C .* (.2969.*sqrt(x./C) - .1260.*(x./C) - ...
    .3516.*(x./C).^2 + .2843.*(x./C).^3 - .1036.*(x./C).^4);
    cp_upper = fnval(Cp_upper,x/C);
    cp_lower = fnval(Cp_lower,x/C);
    p_upper = cp_upper*q + p_inf;
    p_lower = cp_lower*q + p_inf;
    
    N_prime = 0;
    A_prime = 0;

    for j = 1:N 
        N_prime = (N_prime + (-((p_upper(j + 1) + p_upper(j))/2)*((x(j + 1) - x(j)))) + ((p_lower(j + 1) + p_lower(j))/2)*(x(j + 1) - x(j)));
        A_prime = (A_prime + (((p_upper(j + 1) + p_upper(j))/2)*(yt(j + 1) - yt(j))) + ((p_lower(j + 1) + p_lower(j))/2)*(yt(j + 1) - yt(j)));
    end
    Lift = N_prime*cosd(alpha) - A_prime*sind(alpha);
    Drag = N_prime*sind(alpha) + A_prime*cosd(alpha);
end
fprintf('The number of points required for 5 percent relative error is %d\n',N+1);
while ((Lift <= 0) || (((accepted_lift - Lift)/accepted_lift) > .01) && ((accepted_lift - Lift)/accepted_lift) > 0)
    N = N + 1;
    x = linspace(0,C,N + 1);
    yt = (.12/.2) .* C .* (.2969.*sqrt(x./C) - .1260.*(x./C) - ...
    .3516.*(x./C).^2 + .2843.*(x./C).^3 - .1036.*(x./C).^4);
    cp_upper = fnval(Cp_upper,x/C);
    cp_lower = fnval(Cp_lower,x/C);
    p_upper = cp_upper*q + p_inf;
    p_lower = cp_lower*q + p_inf;
    
    N_prime = 0;
    A_prime = 0;

    for j = 1:N 
        N_prime = (N_prime + (-((p_upper(j + 1) + p_upper(j))/2)*((x(j + 1) - x(j)))) + ((p_lower(j + 1) + p_lower(j))/2)*(x(j + 1) - x(j)));
        A_prime = (A_prime + (((p_upper(j + 1) + p_upper(j))/2)*(yt(j + 1) - yt(j))) + ((p_lower(j + 1) + p_lower(j))/2)*(yt(j + 1) - yt(j)));
    end
    Lift = N_prime*cosd(alpha) - A_prime*sind(alpha);
    Drag = N_prime*sind(alpha) + A_prime*cosd(alpha);
end
fprintf('The number of points required for 1 percent relative error is %d\n',N+1);
while ((Lift <= 0) || (((accepted_lift - Lift)/accepted_lift) > .001) && ((accepted_lift - Lift)/accepted_lift) > 0)
    N = N + 1;
    x = linspace(0,C,N + 1);
    yt = (.12/.2) .* C .* (.2969.*sqrt(x./C) - .1260.*(x./C) - ...
    .3516.*(x./C).^2 + .2843.*(x./C).^3 - .1036.*(x./C).^4);
    cp_upper = fnval(Cp_upper,x/C);
    cp_lower = fnval(Cp_lower,x/C);
    p_upper = cp_upper*q + p_inf;
    p_lower = cp_lower*q + p_inf;
    
    N_prime = 0;
    A_prime = 0;

    for j = 1:N 
        N_prime = (N_prime + (-((p_upper(j + 1) + p_upper(j))/2)*((x(j + 1) - x(j)))) + ((p_lower(j + 1) + p_lower(j))/2)*(x(j + 1) - x(j)));
        A_prime = (A_prime + (((p_upper(j + 1) + p_upper(j))/2)*(yt(j + 1) - yt(j))) + ((p_lower(j + 1) + p_lower(j))/2)*(yt(j + 1) - yt(j)));
    end
    Lift = N_prime*cosd(alpha) - A_prime*sind(alpha);
    Drag = N_prime*sind(alpha) + A_prime*cosd(alpha);
end
fprintf('The number of points required for 1/10 percent relative error is %d\n',N+1);

%Plots
%plot of airfoil
figure
hold on
axis equal
plot(x,yt)
plot(x,-yt)
ylabel 'y [m]'
xlabel 'cord length [m]'
title 'NACA 0012 Airfoil'
%Plot of -Cp
figure
hold on
plot(x/C,-cp_upper)
plot(x/C,-cp_lower)
xlabel 'Normalized Chord Length X/C'
ylabel '-Cp'
title 'Cp on NACA 0012 Airfoil'
