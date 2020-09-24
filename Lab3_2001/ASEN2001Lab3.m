clear all
close all
clc

% Define some variables
L = 36*0.0254;
width = 4*0.0254;
Tb = (1/32)*0.0254;
Tf = (3/4)*0.0254;
h = (2*Tb+Tf);
c = h/2;

% Calculate sigma fail
exp_data = xlsread('TestData.xlsx');
exp_data = exp_data(:,2:5);
for i = 2:25
    if exp_data(i-1,4) < exp_data(i-1,2)
        exp_data(i-1,:) = exp_data(i,:);
    end
end
% Define symbols
syms x p0

% Set up the load (per unit length)
qx = width * p0 * sqrt(1-(2*x/L)^2);

% Support reaction 
R = -1 * int(qx,x,-L/2,L/2)/2;

% Shear force and bending moment
Vx = R + int(qx,x,-L/2,x); 
Mx = int(Vx,x,-L/2,x);

% Say moment at x=0 for p0=1
M0_scaled = double(subs(subs(Mx,x,0),p0,1));

% Compute p0 and call p0_calculated
Ib = 1/12*width*(h^3-Tf^3);
If = 1/12*width*Tf^3;
sigma_fail = -7.9773e+06;
sigma_allow = sigma_fail/1.5;
Ef = 0.035; 
Eb = 6.8; 


p0_calculated = (sigma_allow*(Ib+((Ef/Eb)*If)))/(M0_scaled*c);

% Plug in p0_calculated in Mx and Vx
Vx = subs(Vx,p0,p0_calculated);
Mx = subs(Mx,p0,p0_calculated);

% Compute w(x) 
tau_fail = 5.0885e+04;
tau_allow = tau_fail/1.5;

Wx_moment = Mx*c/(sigma_allow*((1/12)*(h^3-Tf^3)+(Ef/Eb)*(1/12)*Tf^3));
Wx_shear = (3/2)*(Vx/(tau_allow*Tf));



% Plot Mx or Vx
dx = -L/2:0.0254:L/2;
figure(1)
subplot(1,2,1)
plot(dx/0.0254,subs(Mx,x,dx))
xlabel('Distance (in)')
ylabel('Shear Force (N)')
subplot(1,2,2)
plot(dx/0.0254,subs(Vx,x,dx))
xlabel('Distance (in)')
ylabel('Bending Moment (Nm)')


w1ans = double(subs(Wx_moment,x,dx))/2;
w2ans = double(abs(subs(Wx_shear,x,dx)))/2;

truew = zeros(1,length(w1ans));
for i = 1:length(w1ans)
   if w1ans(i) > w2ans(i)
       truew(i) = w1ans(i);
   else
       truew(i) = w2ans(i);
   end
end
figure(2)
truew = truew/0.0254;
dx=dx/0.0254;
plot(dx,truew)
hold on
plot(dx, -truew)
hold on
b = linspace(-1.45, 1.45);
plot(-18,b)
hold on
plot(18,b)
ylim([-10,10])
xlabel('Distance (in)')
ylabel('Width (in)')


syms x p0


% Say moment at x=0 for p0=1
M0_scaled = double(subs(subs(Mx,x,0),p0,1));


% Plug in p0_calculated in Mx and Vx
px = p0_calculated*sqrt(1-(2*x/L)^2);


%Find centroids for each region
for i=0:3
    centroid(i+1) = int(x*px,i*L/8,(i+1)*L/8) / int(px,i*L/8,(i+1)*L/8);
    
end

figure(3)
xgrid = -L/2:0.0254:L/2;
plot(xgrid/0.0254,subs(px,x,xgrid))
hold on
a = linspace(0,875);
L = L / 0.0254;
plot(0, a,'k--')
hold on
plot(L/8, a,'k--')
hold on
plot(2*L/8, a,'k--')
hold on
plot(3*L/8 , a,'k--')
hold on
plot(4*L/8 , a,'k--')
hold on
plot(-1*L/8 , a,'k--')
hold on
plot(-2*L/8, a,'k--')
hold on
plot(-3*L/8, a,'k--')
hold on
plot(-4*L/8, a, 'k--')
hold on
for i=1:4
    36*double(centroid(i))
    plot(36*double(centroid(i)), a,'r-')
    plot(-36*double(centroid(i)), a,'r-')
    hold on
end
ylim([0 1000])
xlim([-18 18])
xlabel('Distance (in)')
ylabel('Distributed Pressure (Pa)')

