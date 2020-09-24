clear;close all;clc
addpath('C:\Users\sageh\OneDrive\Documents\MATLAB\ASEN5010\Matlab')
%% CC 12,13
q = [0.1 0.2 0.3]';
qtil = [0 -q(3) q(2);q(3) 0 -q(1);-q(2) q(1) 0];
C = (1/(1+q'*q))*((1 - q'*q)*eye(3) + 2*q*q' - 2*qtil);
clear
BN = [0.333333 -.666667 0.666667;0.871795 0.487179 0.0512821;-0.358974 0.564103 0.74359];
zeta = sqrt(trace(BN) + 1);
q = (1/zeta^2)*[BN(2,3)-BN(3,2);BN(3,1)-BN(1,3);BN(1,2)-BN(2,1)];
clear
qfn = [0.1;0.2;0.3];
qbn = [-0.3;0.3;0.1]; %second rot., q'
qnf = -qfn; %first rot., q
qbf = (qnf-qbn+cross(qnf,qbn))/(1+dot(qnf,qbn));
clear
sigma = [0.1;0.2;0.3];
b0 = (1-(sigma'*sigma)^2)/(1+(sigma'*sigma));
phi = acos(b0)*2;

sigmas = (-sigma)/(sigma'*sigma);
%% CC18
C = MRP2C([0.1 0.2 0.3]');

DCM = [.763314 .0946746 -.639053;-.568047 -.372781 -.733728;-.307692 0.923077 -.230796];

sigma = C2MRP(DCM);
%% CC19

c = addMRP([0.1 0.2 0.3]',[-.1 .3 .1]');

d = addMRP([-.5 -.3 -.1]',[.1 .2 .3]');
%% CC1 att. 
bv1 = [.8273;.5541;-.0920];
bv2 = [-.8285;.5522;-.0955];
nv1 = [-.1517;-.9669;.2050];
nv2 = [-.8393;.4494;-.3044];

b1 = bv1;
b2 = cross(bv1,bv2)/norm(cross(bv1,bv2));
b3 = cross(b1,b2);

n1 = nv1;
n2 = cross(nv1,nv2)/norm(cross(nv1,nv2));
n3 = cross(n1,n2);

bDCM = [b1 b2 b3];
nDCM = [n1 n2 n3];

DCM = bDCM*nDCM';

BbarN = [0.969846 0.17101 0.173648;-0.200706 0.96461 0.17101;-0.138258 -0.200706 0.969846]

BN = [0.963592 0.187303 0.190809;
    -.223042 0.956645 0.187303;
    -.147454 -.223042 .963592];

DCM = BN*BbarN';

phi = acosd((.5*(trace(DCM) - 1)))
%% CC 2 - Kinetic Energy
R1 = [1;-1;2];
R2 = [-1;-3;2];
R3 = [2;-1;-1];
R4 = [3;-1;-2];
R1dot = [2;1;1];
R2dot = [0;-1;1];
R3dot = [3;2;-1];
R4dot = [0;0;1];

I1 = [R1(2)^2+R1(3)^2 -R1(1)*R1(2) -R1(1)*R1(3);
        -R1(1)*R1(2) R1(1)^2+R1(3)^2 -R1(2)*R1(3);
        -R1(1)*R1(3) -R1(2)*R1(3) R1(2)^2+R1(2)^2];
I2 = [R2(2)^2+R2(3)^2 -R2(1)*R2(2) -R2(1)*R2(3);
    -R2(1)*R2(2) R2(1)^2+R2(3)^2 -R2(2)*R2(3);
    -R2(1)*R2(3) -R2(2)*R2(3) R2(2)^2+R2(2)^2];
I3 = [R3(2)^2+R3(3)^2 -R3(1)*R3(2) -R3(1)*R3(3);
    -R3(1)*R3(2) R3(1)^2+R3(3)^2 -R3(2)*R3(3);
    -R3(1)*R3(3) -R3(2)*R3(3) R3(2)^2+R3(2)^2];
I4 = [R4(2)^2+R4(3)^2 -R4(1)*R4(2) -R4(1)*R4(3);
    -R4(1)*R4(2) R4(1)^2+R4(3)^2 -R4(2)*R4(3);
    -R4(1)*R4(3) -R4(2)*R4(3) R4(2)^2+R4(2)^2];

m1 = 1;
m2 = 1;
m3 = 2;
m4 = 2;

w1 = R1dot./R1;
w2 = R2dot./R2;
w3 = R3dot./R3;
w4 = R4dot./R4;

T_trans = .5*(m1*dot(R1dot,R1dot)+m2*dot(R2dot,R2dot)+m3*dot(R3dot,R3dot)+m4*dot(R4dot,R4dot));

T_rot = .5*(w1'*I1*w1 + w2'*I2*w2 + w3'*I3*w3 + w4'*I4*w4);

COM = (R1dot*1 + R2dot*1 + R3dot*2 + R4*2)/6;

%% CC5
I = [10 1 -1;1 5 1;-1 1 8];
w = [0.01;-0.01;0.01];
C = Euler3212C([deg2rad(-10),deg2rad(10),deg2rad(5)]');

