clear all;close all;clc
%% Question 3
% part a
AU = 1.5e8; %km
mu_s = 1.33e11; %mu sun
mu_e = 4e5; %mu earth
Re = 6400; %radius of earth
aorig = 1.5*AU;
eorig = 0.5;
%calculate Tisserand's Crit. prior to Earth encounter
n = sqrt(mu_s/(AU^3));
Jconst = -mu_s/(2*aorig) - n*sqrt(mu_s*aorig*(1-eorig^2));
% syms a
% e = [0:.02:1];
% for i = 1:numel(e)
%     J = -mu_s/(2*a) - n*sqrt(mu_s*a*(1-e(i)^2));
%     soln(i) = vpasolve(J==Jconst);
% end
a = linspace(6e7,5e9,200);
% e = sqrt(1 - ((-(Jconst/mu_s) - (1/2.*a)).*(AU^(3/2))).^2.*(1./a));
e = sqrt(1-((-Jconst/mu_s-1./(2*a))*(AU^(3/2))).^2./a);
figure
hold on
plot(a,e)
ylim([0 1])
ylabel('Eccentricity')
xlabel('Semi-Major Axis')
title('a VS e')
grid on
grid minor
% part b
V_SC0 = sqrt(mu_s*((2/AU) - 1/(1.5*AU)));
Vp = sqrt(mu_s/AU);
f = acos((((1.5*AU)*(1 - .5^2))/AU - 1)*(1/.5));
gamma = atan((.5*sin(f))/(1 + .5*cos(f)));
v_inf = sqrt(Vp^2 + V_SC0^2 - 2*Vp*V_SC0*cos(gamma));

theta = asin((V_SC0*sin(gamma))/v_inf);
delta_star = 2*asin(1/(sqrt(1 + ((v_inf^2*Re)/mu_e))));
thetaprimep = theta + delta_star;
thetaprimen = theta - delta_star;

V_SCfp = sqrt(Vp^2 + v_inf^2 - 2*Vp*v_inf*cos(thetaprimep));
V_SCfn = sqrt(Vp^2 + v_inf^2 - 2*Vp*v_inf*cos(thetaprimen));

aprime_upper = 1/((2/AU) - (V_SCfp^2/mu_s));
aprime_lower = 1/((2/AU) - (V_SCfn^2/mu_s));

elow = sqrt(1-((-Jconst/mu_s-1./(2*aprime_lower))*(AU^(3/2))).^2./aprime_lower);
ehigh = sqrt(1-((-Jconst/mu_s-1./(2*aprime_upper))*(AU^(3/2))).^2./aprime_upper);
plot(aorig,eorig,'*')
plot(aprime_lower,elow,'*')
plot(aprime_upper,ehigh,'*')
legend('Tisserand Criterion','Original Orbit Parameters','Lower Bound of Orbit Parameters','Upper Bound of Orbit Paramters')

V_SCf = sqrt(Vp^2 + v_inf^2 - 2*Vp*v_inf*cos(delta_star));
dv = V_SCf - V_SC0;
inc = 2*asin((-dv/(2*V_SC0)));
% inc = deg2rad(5);
%recalc a' and e based on new inclination
aprime = 1/((2/AU) - (V_SCf^2/mu_s));
e_inc = sqrt(1-((-Jconst/mu_s-1./(2*aprime))*(AU^(3/2))/cos(0.0873)).^2./aprime);
plot(aprime,e_inc,'*')
legend('Tisserand Criterion','Original Orbit Paramters','Lower Bound of Orbit Parameters','Upper Bound of Orbit Paramters','Parameters after Max Inclination')
%% Question 5, part b
x = linspace(-2,2,200);
y = linspace(-2,2,200);
z = 0;

for i = 1:numel(x)
    for j = 1:numel(y)
        Vh = .5*(3*x(i)^2 - z^2) + (1/(sqrt(x(i)^2 + y(j)^2 + z^2)));
        J(i,j) = -Vh;
    end
end
figure
contour(x,y,J,200)
xlabel('X')
ylabel('Y')
title('Constant Jacobi Constant Values, X-Y Plane')
hold on
plot(1/(3^(1/3)),0,'*')
plot(-1/(3^(1/3)),0,'*')

x = linspace(-2,2,100);
y = 0;
z = linspace(-2,2,100);
for i = 1:numel(x)
    for j = 1:numel(z)
        Vh = .5*(3*x(i)^2 - z(j)^2) + (1/(sqrt(x(i)^2 + y^2 + z(j)^2)));
        J2(i,j) = -Vh;
    end
end
figure
contour(x,z,J2,200)
xlabel('X')
ylabel('Z')
title('Constant Jacobi Constant Values, X-Z Plane')
hold on
plot(1/(3^(1/3)),0,'*')
plot(-1/(3^(1/3)),0,'*')