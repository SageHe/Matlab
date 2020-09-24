clear all;close all;clc
%input given three position vectors
r1 = [-1.054e8 1.579e8 -1.520e5];
r2 = [-1.461e8 1.081e8 -2.265e5];
r3 = [-1.652e8 4.254e7 -2.673e5];

magr1 = norm(r1);
magr2 = norm(r2);
magr3 = norm(r3);
%calculate variables C12, C23, and C31
C12 = cross(r1,r2);
C23 = cross(r2,r3);
C31 = cross(r3,r1);
%Calculate N,D,and S
N = magr1.*C23 + magr2.*C31 + magr3.*C12;
D = C12 + C23 + C31;
S = r1.*(magr2 - magr3) + r2.*(magr3 - magr1) + r3.*(magr1 - magr2); 
%determine v2
v2 = sqrt((132712000000)./(norm(N).*norm(D))).*(((cross(D,r2))./magr2) + S);
vr = (dot(v2,r2)/norm(r2));
e = (1/132712000000)*((norm(v2)^2 - (132712000000/norm(r2))).*r2 - norm(r2)*vr.*v2);