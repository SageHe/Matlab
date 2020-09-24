%Sage Herrin
%ASEN 3128 HW8
%SID 106071909
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Homework 7 code to be used 
clear all;close all;clc
%define constants
global ss theta_0 u_0 A mat
g = 9.81;
S = 510.9667;%m^2
b = 59.6433;%m
c_bar = 8.3241;%m
theta_0 = 0;
W = 2.8317e6;%Newtons
m = W/9.81;
V = 157.886;%m/s
rho = .6532;%kg/m^3
u_0 = V;%m/s
xi = deg2rad(-6.8);
C_w0 = W/(.5*rho*(u_0^2)*S);%weight coeff
I_x = 24675586.69;%kgm^2
I_y = 44877574.145;%kgm^2
I_z = 67384152.115;%kgm^2

%non-dimensional stability matrix
ndstab = [-0.1080 -0.1060 0.1043;0.2193 -4.920 -1.023;0 -5.921 -23.92;0 5.896 -6.314];
%Create dimensional stability matrix for body frame using above constants 
bs = [rho*u_0*S*C_w0*sin(theta_0) + .5*rho*u_0*S*ndstab(1,1) -rho*u_0*S*C_w0*cos(theta_0) + .5*rho*u_0*S*ndstab(1,2) .5*rho*u_0*c_bar*S*ndstab(1,3);...
        .5*rho*u_0*S*ndstab(2,1) .5*rho*u_0*S*ndstab(2,2) .5*rho*u_0*c_bar*S*ndstab(2,3);...
        .25*rho*u_0*c_bar*S*ndstab(3,1) .25*rho*u_0*c_bar*S*ndstab(3,2) .25*rho*u_0*(c_bar^2)*S*ndstab(3,3);...
        .25*rho*c_bar*S*ndstab(4,1) .25*rho*c_bar*S*ndstab(4,2) .25*rho*(c_bar^2)*S*ndstab(4,3)];
    
%change above matrix from body frame to stability frame using following
%conversion, taken from Appendix B12 of textbook
ss = zeros(4,3);
ss(1,1) = bs(1,1)*cos(xi)^2 - (bs(2,1) + bs(1,2))*sin(xi)*cos(xi) + bs(2,2)*sin(xi)^2;
ss(2,1) = bs(2,1)*cos(xi)^2 + (bs(1,1) - bs(2,2))*sin(xi)*cos(xi) - bs(1,2)*sin(xi)^2;
ss(3,1) = bs(3,1)*cos(xi) - bs(3,2)*sin(xi);
ss(4,1) = -bs(4,2)*sin(xi)*cos(xi);
ss(1,2) = bs(1,2)*cos(xi)^2 - (bs(2,2) - bs(1,1))*sin(xi)*cos(xi) - bs(2,1)*sin(xi)^2;
ss(2,2) = bs(2,2)*cos(xi)^2 + (bs(1,2) + bs(2,1))*sin(xi)*cos(xi) + bs(1,1)*sin(xi)^2;
ss(3,2) = bs(3,2)*cos(xi) + bs(3,1)*sin(xi);
ss(4,2) = bs(4,2)*cos(xi)^2;
ss(1,3) = bs(1,3)*cos(xi) - bs(2,3)*sin(xi);
ss(2,3) = bs(2,3)*cos(xi) + bs(1,3)*sin(xi);
ss(3,3) = bs(3,3);
ss(4,3) = bs(4,3)*cos(xi);
%the above stability matrix converted to the stability frame is used to
%construct the A matrix of the linearized longitudinal model
A = [(ss(1,1)/m) (ss(2,1)/m) 0 -g*cos(theta_0);...
    (ss(1,2)/(m - ss(4,2))) (ss(2,2)/(m - ss(4,2))) ((ss(3,2) + m*u_0)/(m - ss(4,2))) ((-m*g*sin(theta_0))/(m - ss(4,2)));...
    ((1/I_y)*(ss(1,3) + ((ss(4,3)*ss(1,2))/(m - ss(4,2))))) (1/I_y)*(ss(2,3) + ((ss(4,3)*ss(2,2))/(m - ss(4,2))))...
    ((1/I_y)*(ss(3,3) + ((ss(4,3)*(ss(3,2) + m*u_0)/(m - ss(4,2)))))) (-(ss(4,3)*m*g*sin(theta_0))/(I_y*(m - ss(4,2))));...
    0 0 1 0];
    
%determine the eigen values and eigenvectors of A using built in eig
%function
[V,D] = eig(A);

%Calculate the natural frequencies and damping ratios using built in
%function damp
[Wn,Z] = damp(D);

%compare the eigenvalues determined above to the approx given in class and
%compare oscillation period of phugoid mode above to Lanchester approx. in
%textbook
%given approx of short period mode evalues
eval_SP_1 = (ss(3,3)/(2*I_y)) + sqrt((ss(3,3))^2 + (4*I_y*u_0*(ss(4,3))))*(1/(2*I_y));
eval_SP_2 = (ss(3,3)/(2*I_y)) - sqrt((ss(3,3))^2 + (4*I_y*u_0*(ss(4,3))))*(1/(2*I_y));

%calculated phugoid period from previous values
T_ph_1 = (2*pi)/0.0938;
%phuoid oscillation period from textbook
T_ph_2 = 0.138*u_0;

%simulate linearized longitudinal dynamics to verify trim state is an
%equilibrium and perturb initial states one at a time
delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 500];
[t,y] = ode45('hw7ode',tspan,inish_condish);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 3a, design for increased pitch stiffness in short period mode by a
%variable scale factor Ks ranging from 1-3, in increments of 0.01, while
%retaining orig. damping ratio, and use short period approx, from class
%define constants to be used 
X_delta_e = -15.883;
Z_delta_e = -1.5176e6;
M_delta_e = -5.00035e7;
zeta = .3862;
%solve for K2 using pitch stiffness
Ks = [1:.01:3];
for i = 1:length(Ks)
    K2 = (u_0*ss(2,3)*(1 - Ks(i)))/(M_delta_e);
    Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
    K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
    %calculate new A matrix with controls and plot evalues
    A = [(ss(3,3)/I_y) - (M_delta_e*K1)/I_y ((u_0*ss(2,3))/I_y) - ((M_delta_e*K2)/I_y); 1 0];
    [V,D] = eig(A);
    eig1(i) = D(1,1);
    eig2(i) = D(2,2);
end
%calculate new A matrix with controls and plot evalues

figure(1)
hold on
grid on
grid minor
plot(eig1,'o')
plot(eig2,'o')
title('Locus of Eigenvalues')
legend('\lambda 1','\lambda 2')

%Part b, implement range of controls laws on full linearized longitudinal
%dynamics and calc. actual short period and phugoid mode evalues, plot in
%complex plane.

A = [(ss(1,1)/m) (ss(2,1)/m) 0 -g*cos(theta_0);...
    (ss(1,2)/(m - ss(4,2))) (ss(2,2)/(m - ss(4,2))) ((ss(3,2) + m*u_0)/(m - ss(4,2))) ((-m*g*sin(theta_0))/(m - ss(4,2)));...
    ((1/I_y)*(ss(1,3) + ((ss(4,3)*ss(1,2))/(m - ss(4,2))))) (1/I_y)*(ss(2,3) + ((ss(4,3)*ss(2,2))/(m - ss(4,2))))...
    ((1/I_y)*(ss(3,3) + ((ss(4,3)*(ss(3,2) + m*u_0)/(m - ss(4,2)))))) (-(ss(4,3)*m*g*sin(theta_0))/(I_y*(m - ss(4,2))));...
    0 0 1 0];

B = [-5.5025e-5 0;-5.3322 0;-1.1142 0;0 0];

Ks = [1:.01:3];
for i = 1:length(Ks)
    K2 = (u_0*ss(2,3)*(1 - Ks(i)))/(M_delta_e);
    Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
    K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
    %calculate new A matrix with controls and plot evalues
    K = [0 0 -K1 -K2;0 0 0 0];
    mat = A + B*K;
    [V,D] = eig(mat);
    %sort short period into eig1 and eig2, phugoid at eig3 and eig4
    if abs(real(D(1,1))) > abs(real(D(3,3)))
        eig1(i) = D(1,1);
        eig2(i) = D(2,2);
        eig3(i) = D(3,3);
        eig4(i) = D(4,4);
        if Ks(i) == 1
            figure(2)
            hold on
            grid on
            grid minor
            plot(eig1(i),'*','MarkerSize',20)
            plot(eig2(i),'*','MarkerSize',20)
%             legend('Ks=1','Ks=1')
            figure(3)
            hold on
            grid on
            grid minor
            plot(eig3(i),'*','MarkerSize',20)
            plot(eig4(i),'*','MarkerSize',20)
%             legend('Ks=1','Ks=1')
        end
        if Ks(i) == 2
            figure(2)
            hold on
            grid on
            grid minor
            plot(real(eig1(i)),imag(eig1(i)),'*','MarkerSize',20)
            plot(real(eig2(i)),imag(eig2(i)),'*','MarkerSize',20)
            legend('Ks=1','Ks=1','Ks=2','Ks=2')
            figure(3)
            hold on
            grid on
            grid minor
            plot(real(eig3(i)),imag(eig3(i)),'*','MarkerSize',20)
            plot(real(eig4(i)),imag(eig4(i)),'*','MarkerSize',20)
            legend('Ks=1','Ks=1','Ks=2','Ks=2')
        end
    else
        eig1(i) = D(3,3);
        eig2(i) = D(4,4);
        eig3(i) = D(1,1);
        eig4(i) = D(2,2);
        if Ks(i) == 1
            figure(2)
            hold on
            grid on
            grid minor
            plot(eig1(i),'*')
            plot(eig2(i),'*')
            legend('Ks=1','Ks=1')
            figure(3)
            hold on
            grid on
            grid minor
            plot(eig3(i),'*')
            plot(eig4(i),'*')
            legend('Ks=1','Ks=1')
        end
        if Ks(i) == 2
            figure(2)
            hold on
            grid on
            grid minor
            plot(real(eig1(i)),imag(eig1(i)),'*')
            plot(real(eig2(i)),imag(eig2(i)),'*')
            legend('Ks=2','Ks=2')
            figure(3)
            hold on
            grid on
            grid minor
            plot(real(eig3(i)),imag(eig3(i)),'*')
            plot(real(eig4(i)),imag(eig4(i)),'*')
            legend('Ks=2','Ks=2')
        end
    end
end

figure(2)
hold on
grid on
grid minor
plot(eig1,'o')
plot(eig2,'o')
title('Short Period Mode')
xlabel('Real')
ylabel('Imaginary')
figure(3) 
hold on
grid on
grid minor
plot(eig3,'o')
plot(eig4,'o')
title('Phugoid Mode')
xlabel('Real')
ylabel('Imaginary')

Ks = 1;
K2 = (u_0*ss(2,3)*(1 - Ks))/(M_delta_e);
Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
%calculate new A matrix with controls and plot evalues
K = [0 0 -K1 -K2;0 0 0 0];
mat = A + B*K;
[V,D] = eig(mat);
% figure(4)
% hold on
% plot(D(1,1),'*')
% plot(D(2,2),'*')
% title('short period')
% legend('Ks1','Ks1')
% figure(5)
% hold on
% plot(D(3,3),'*')
% plot(D(4,4),'*')
% legend('Ks1','Ks1')

delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 200];
[t,y] = ode45('hw8ode',tspan,inish_condish);

figure(4)
suptitle('Modal Response When Ks=1')
subplot(2,3,1)
plot(t,y(:,1))
grid on
grid minor
title('\Delta u')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,2)
plot(t,y(:,2))
grid on
grid minor
title('\Delta w')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,3)
plot(t,y(:,3))
grid on 
grid minor
title('\Delta q')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,4)
plot(t,y(:,4))
grid on
grid minor
title('\Delta \theta')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,5)
plot(t,y(:,5))
grid on
grid minor
title('\Delta x')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,6)
plot(t,y(:,6))
grid on
grid minor
title('\Delta z')
xlabel('Time (s)')
ylabel('Amplitude')

Ks = 2;
K2 = (u_0*ss(2,3)*(1 - Ks))/(M_delta_e);
Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
%calculate new A matrix with controls and plot evalues
K = [0 0 -K1 -K2;0 0 0 0];
mat = A + B*K;
[V,D] = eig(mat);
% figure(4)
% hold on
% plot(D(1,1),'*')
% plot(D(2,2),'*')
% legend('Ks2','Ks2')
% figure(5)
% hold on
% plot(D(3,3),'*')
% plot(D(4,4),'*')
% legend('Ks2','Ks2')

delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 200];
[t,y] = ode45('hw8ode',tspan,inish_condish);

figure(5)
suptitle('Modal Response When Ks=2')
subplot(2,3,1)
plot(t,y(:,1))
grid on
grid minor
title('\Delta u')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,2)
plot(t,y(:,2))
grid on
grid minor
title('\Delta w')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,3)
plot(t,y(:,3))
grid on
grid minor
title('\Delta q')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,4)
plot(t,y(:,4))
grid on
grid minor
title('\Delta \theta')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,5)
plot(t,y(:,5))
grid on
grid minor
title('\Delta x')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,3,6)
plot(t,y(:,6))
grid on
grid minor
title('Delta z')
xlabel('Time (s)')
ylabel('Amplitude')
%question d, original modal response

Ks = 2;
K2 = (u_0*ss(2,3)*(1 - Ks))/(M_delta_e);
Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
%calculate new A matrix with controls and plot evalues
K = [0 0 -K1 -K2;0 0 0 0];
mat = A + B*K;

delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 200];
[t1,y1] = ode45('hw8ode',tspan,inish_condish);

%modal response with decoupling in control matrix B
B1 = B;
B1([1 2]) = 0;

Ks = 2;
K2 = (u_0*ss(2,3)*(1 - Ks))/(M_delta_e);
Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
%calculate new A matrix with controls and plot evalues
K = [0 0 -K1 -K2;0 0 0 0];
mat = A + B1*K;

delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 200];
[t2,y2] = ode45('hw8ode',tspan,inish_condish);

%modal response with decoupling in A matrix 
A1 = A;
A1([1;2],[1:2]) = 0;
A1([3:4],[3:4]) = 0;

Ks = 2;
K2 = (u_0*ss(2,3)*(1 - Ks))/(M_delta_e);
Wn = sqrt(((-u_0*ss(2,3))/(I_y)) + ((M_delta_e*K2)/I_y));
K1 = ((2*zeta*Wn + (ss(3,3)/(I_y)))*(I_y/M_delta_e));
%calculate new A matrix with controls and plot evalues
K = [0 0 -K1 -K2;0 0 0 0];
mat = A1 + B*K;

delta_u = 0; %m/s
delta_w = 0; %m/s
delta_q = 0; %rad/s
delta_theta = 0.1; %rad/s
delta_x = 0;
delta_z = 0;
inish_condish = [delta_u delta_w delta_q delta_theta delta_x delta_z];
tspan = [0 200];
[t3,y3] = ode45('hw8ode',tspan,inish_condish);

figure(6)
suptitle('Comparison Between Decoupling in A and B Matrix')
subplot(2,3,1)
hold on
grid on
grid minor
plot(t1,y1(:,1))
plot(t2,y2(:,1))
plot(t3,y3(:,1))
title('\Delta u')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')
subplot(2,3,2)
hold on
grid on
grid minor
plot(t1,y1(:,2))
plot(t2,y2(:,2))
plot(t3,y3(:,2))
title('\Delta w')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')
subplot(2,3,3)
hold on
grid on
grid minor
plot(t1,y1(:,3))
plot(t2,y2(:,3))
plot(t3,y3(:,3))
title('\Delta q')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')
subplot(2,3,4)
hold on
grid on
grid minor
plot(t1,y1(:,4))
plot(t2,y2(:,4))
plot(t3,y3(:,4))
title('\Delta \theta')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')
subplot(2,3,5)
hold on
grid on
grid minor
plot(t1,y1(:,5))
plot(t2,y2(:,5))
plot(t3,y3(:,5))
title('\Delta x')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')
subplot(2,3,6)
hold on
grid on
grid minor
plot(t1,y1(:,6))
plot(t2,y2(:,6))
plot(t3,y3(:,6))
title('\Delta z')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Orig','B decouple','A decouple')


