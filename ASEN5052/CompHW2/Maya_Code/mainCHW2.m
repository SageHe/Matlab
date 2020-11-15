%% ASEN 5519 Comp Hw 2/Analytical HW 4
clear; close all;

%% 1ci
J2 = 1e-3;
RE = 6400;
mu = 4e5;
O_dot = 2*pi/24/365/3600;
C = O_dot*2/3/J2/RE^2/sqrt(mu);
aM = C^(-2/7);

%% 1cii
xa = linspace(6400,aM,100);
yi = acos(-(xa/aM).^(7/2));

figure(1)
plot(xa,yi*180/pi)
    xlabel('Orbit Radius (km)');
    ylabel('Inclination (deg)');
    title('Inlincation for a Sun-Synchronous Orbit vs Orbit Radius');
grid on

%% 1ciii
Tmax = 2*pi*sqrt(aM^3/mu);
Tmin = 2*pi*sqrt(RE^3/mu);

%% 1di

rp = 7000;
i = asin(2/sqrt(5));
T = 12*3600;
a = (T^2*mu/4/pi^2)^(1/3);
e = (2*a - 2*rp)/(2*a);
n = sqrt(mu/(a^3));
rate = -3*J2*RE^2*n*cos(i)/2/a^2/(1-e^2)^2;


rate_deg = rate*180/pi;
time_1deg = 1/rate_deg;

%% 1dii
f = pi/2;
E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
tHalfSouth = sqrt(a^3/mu)*(E - e*sin(E));
tNorth = T - 2*tHalfSouth;
percentNorth = tNorth/T;


%% 1diii

i_prime = i - .5*pi/180;
womega_dot = 3*J2*RE^2*n*(2 - (5*sin(i_prime)^2/2))/2/a^2/(1-e^2)^2;
t180 = pi/womega_dot;


%% 1ei
i = 0;
z = 200;
a = RE + z;
n = sqrt(mu/a^3);
nBar = n*(1 + 3*J2*RE^2/2/a^2);
error = (nBar-n)/nBar;

tError180 = pi/(nBar-n);

%% 1eii

d = (nBar - n)*2*pi*sqrt(a^3/mu);
arc = d*a;



%% 2a
% using ODE45
a = 7000;
e = 0.1;
i = 45*pi/180;
womega = 0;
Omega = 0;
t = 0;
mu = 4e5;
[rvec0, vvec0] = rvElements(a,e,i,womega,Omega,t,mu);

x = rvec0(1);
y = rvec0(2);
z = rvec0(3);
u = vvec0(1);
v = vvec0(2);
w = vvec0(3);
mu = 4e5;

T = 2*pi*sqrt(a^3/mu);
IC = [x y z u v w];
var = [mu, RE, J2];
% gives pos and velocity over specified interval (30 orbits for this)
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,r] = ode45(@(t,IC) analyticalJ2(t,IC,var),[0 100*T],IC,opts);

% plotting
figure(2)
plot3(r(:,1),r(:,2),r(:,3))
    xlabel('X Position (km)');
    ylabel('Y Position (km)');
    zlabel('Z Position (km)');
    title('Position Over 100 Orbits With J2 Perturbation');
grid on



% H
hvecs = cross(r(:,1:3),r(:,4:6));
hvecs = hvecs(:,3);
figure(3)
plot(t,hvecs)
    xlabel('Time (s)');
    ylabel('Specific Angular Momentum (km^2/s)');
    title('Specific Angular Momentum Over Time');
    ylim([hvecs(1)-1, hvecs(end)+1])
grid on

% E
Evecs = zeros(length(hvecs),1);
for j = 1:length(hvecs)
	v = norm(r(j,4:6));
    radi = norm(r(j,1:3));
    z = r(j,3);
    Evecs(j) = .5*v^2 - mu/radi - mu*RE^2*J2*(1-3*z^2/radi^2)/2/radi^3;
end

figure(4)
plot(t,Evecs)
    xlabel('Time (s)');
    ylabel('Specific Energy (km^2/s^2)');
    title('Specific Energy Over Time');
    ylim([Evecs(1)-1, Evecs(end)+1])
grid on


%% 2b

[t,r] = ode45(@(t,IC) analyticalJ2(t,IC,var),[0 10*T],IC,opts);
for j = 1:length(t)
    % numerical
    [~,~,~,elems] = rvPredict(r(j,1:3),r(j,4:6),t(j),mu);
    anum(j) = elems(1);
    enum(j) = elems(2);
    inum(j) = elems(3)*180/pi;
    wnum(j) = elems(4)*180/pi;
    Onum(j) = elems(5)*180/pi;
    signum(j) = elems(6) - enum(j)*sin(elems(6));
    
    % averaging theory
    aavg(j) = a;
    eavg(j) = e;
    iavg(j) = i*180/pi;
    n = sqrt(mu/a^3);
    p = a*(1-e^2);
    wavg(j) = womega*180/pi + 3*J2*RE^2*n*(2-5*sin(i)^2/2)*t(j)/2/p^2*180/pi/2;
    Oavg(j) = Omega*180/pi - 3*J2*RE^2*n*cos(i)*t(j)/2/p^2*180/pi/2;
    sigavg(j) = 3*J2*RE^2*n*(1-3*sin(i)^2/2)*sqrt(1-e^2)/2/p^2;
end
    
figure(5)
subplot(2,3,1)
plot(t,anum)
hold on
plot(t,aavg)
        xlabel('Time (s)')
        ylabel('(km)')
        title('Semi-Major Axis, a')
        legend('numerical','averaged')
        ylim([anum(1)-10 anum(1)+10])

subplot(2,3,2)
plot(t,enum)
hold on
plot(t,eavg)
        xlabel('Time (s)')
        legend('numerical','averaged')
        ylim([enum(1)-.01 enum(1)+.01])
        title('Eccentricity, e')

subplot(2,3,3)
plot(t,inum)
hold on
plot(t,iavg)
ylim([iavg(1)-.1 iavg(1)+.1])
    xlabel('Time (s)')
    ylabel('(deg)')
    legend('numerical','averaged')
    title('Inclination, i')


subplot(2,3,4)
plot(t,wnum)
hold on
plot(t,wavg)
    ylabel('(deg)')
    xlabel('Time (s)')
    title('Argument of Periapsis, \omega')
    % ylim([wavg(1)-.1 wavg(end)+.1])
    legend('numerical','averaged')


subplot(2,3,5)
plot(t,Onum)
hold on
plot(t,Oavg)
    xlabel('Time (s)')
    ylabel('(deg)')
    title('Right Ascension of Ascending Node, \Omega')
    % ylim([Oavg(1)-.1 Oavg(end)+.1])
    legend('numerical','averaged')


subplot(2,3,6)
plot(t,signum)
hold on
plot(t,sigavg)
    xlabel('Time (s)')
    ylabel('(deg)')
    title('\sigma')
    legend('numerical','averaged')
    
set(gcf, 'Position',  [100, 100, 950, 500])

%% 4a

mu = 1;
a = 1;
i = 0;
e = 0;
Omega = 0;
womega = 0;
t = 0;

[rvec4, vvec4] = rvElements(a,e,i,womega,Omega,t,mu);

g = [.001 .01 .1 1];
IC = [rvec4 vvec4];

T = 2*pi*sqrt(1/mu);
for j = 1:4
    U = [g(j) 0 0];
    var = {mu,RE,U};
    [time{j},sol{j}] = ode45(@(t,IC) analyticalg(t,IC,var),[0 30*T],IC,opts);
end

% h about x
for j = 1:4
    h = cross(sol{j}(:,1:3),sol{j}(:,4:6));
    hx{j} = h(:,1);

figure(5+j)
plot(time{j},hx{j})
    xlabel('Time (s)');
    ylabel('Specific Angular Momentum About the x-axis (km^2/s)');
    titleString = sprintf('Specific Angular Momentum Over Time, g = %1.3f',g(j));
    title(titleString);
grid on

end


% energy
for j = 1:4
    for k = 1:length(sol{j}(:,1))
        radi = norm(sol{j}(k,1:3));
        v = norm(sol{j}(k,4:6));
        x = sol{j}(k,1);
        Energy{j}(k) = .5*v^2 - mu/radi - g(j)*x;
    end
    figure(9+j)
    plot(time{j},Energy{j})
    xlabel('Time (s)');
    ylabel('Specific Energy (km^2/s^2)');
    titleString = sprintf('Specific Energy Over Time, g = %1.3f',g(j));
    title(titleString);
    ylim([Energy{j}(1)-1, Energy{j}(1)+1])
    grid on
end

%% 4b

% plot trajectory
for j = 1:4
    
    figure(13+j)
    plot3(sol{j}(:,1),sol{j}(:,2),sol{j}(:,3))
        xlabel('X Posisition (km)');
        ylabel('Y Position (km)');
        zlabel('Z Position (km)');
        titleString = sprintf('Trajectory Over 30 Orbits, g = %1.3f',g(j));
        title(titleString);
    grid on
    view([.6,.4,1])

end

% a,e,wt
clear anum enum wnum Onum aavg eavg
for j = 1:4
    for k = 1:length(time{j})
    % numerical
    [~,~,~,elems] = rvPredict(sol{j}(k,1:3),sol{j}(k,4:6),time{j}(k),mu);
    anum{j}(k) = elems(1);
    enum{j}(k) = elems(2);
    wnum = elems(4)*180/pi;
    Onum = elems(5)*180/pi;
    wt{j}(k) = wnum+Onum;
    
    % averaging theory
    aavg{j}(k) = a;
    eavg{j}(k)= 0;
    wtavg{j}(k) = 0;
    end
    figure(17+j)
    subplot(3,1,1)
        plot(time{j},anum{j})
        hold on
        plot(time{j},aavg{j})
        xlabel('Time (s)')
        ylabel('(km)')
        title('Semi-Major Axis, a')
        legend('numerical','averaged')
        ylim([anum{j}(1)-10 anum{j}(1)+10])
    subplot(3,1,2)
        plot(time{j},enum{j})
        hold on
        plot(time{j},eavg{j})
        xlabel('Time (s)')
        title('Eccentricity, e')
        legend('numerical','averaged')
%         ylim([anum{j}(1)-10 anum{j}(1)+10])
    subplot(3,1,3)
        plot(time{j},wt{j})
        hold on
        plot(time{j},wtavg{j})
        xlabel('Time (s)')
        ylabel('(deg)')
        title('Longitude of Periapsis, \omega ~')
        legend('numerical','averaged')
%         ylim([anum{j}(1)-10 anum{j}(1)+10])
end



%% save figures
% saveas(figure(1),'sunSynci.png')
% saveas(figure(2),'j2pert.png')
% saveas(figure(3),'j2angmom.png')
% saveas(figure(4),'j2energy.png')
% saveas(figure(5),'j2numavg.png')
% saveas(figure(6),'angmomg001.png')
% saveas(figure(7),'angmomg01.png')
% saveas(figure(8),'angmomg1.png')
% saveas(figure(9),'angmomg10.png')
% saveas(figure(10),'energyg001.png')
% saveas(figure(11),'energyg01.png')
% saveas(figure(12),'energyg1.png')
% saveas(figure(13),'energyg10.png')
% saveas(figure(14),'posg001.png')
% saveas(figure(15),'posg01.png')
% saveas(figure(16),'posg1.png')
% saveas(figure(17),'posg10.png')
% saveas(figure(18),'elemg001.png')
% saveas(figure(19),'elemg01.png')
% saveas(figure(20),'elemg1.png')
% saveas(figure(21),'elemg10.png')





















