%% 3.17 
clear all;close all;clc
phi = 45;
e = [1 0 0]';

etil = [0 -e(3) e(2);e(3) 0 -e(1);-e(2) e(1) 0];

C = expm(-deg2rad(phi)*etil);

compC = eye(3)*cosd(phi) - sind(phi)*etil + (1 - cosd(phi))*e*transpose(e);
%% 3.22
clear all;close all;clc
syms e1 e2 e3 phi b0 b1 b2 b3
c = cos(phi);
s = sin(phi);
sig = 1 - c;

C = [e1^2*sig+c e1*e2*sig+e3*s e1*e3*sig-e2*s;e2*e1*sig-e3*s e2^2+c e2*e3*sig+e1*s;e3*e1*sig+e2*s e3*e2*sig-e1*s e3^2+c];

c = 2*cos(phi/2)^2 - 1;
s = 2*sin(phi/2)*cos(phi/2);
sig = 1 - c;

C = [e1^2*sig+c e1*e2*sig+e3*s e1*e3*sig-e2*s;e2*e1*sig-e3*s e2^2*sig+c e2*e3*sig+e1*s;e3*e1*sig+e2*s e3*e2*sig-e1*s e3^2*sig+c];

% C =
%  
% [            2*cos(phi/2)^2 - e1^2*(2*cos(phi/2)^2 - 2) - 1,   2*e3*cos(phi/2)*sin(phi/2) - e1*e2*(2*cos(phi/2)^2 - 2), - 2*e2*cos(phi/2)*sin(phi/2) - e1*e3*(2*cos(phi/2)^2 - 2)]
% [ - 2*e3*cos(phi/2)*sin(phi/2) - e1*e2*(2*cos(phi/2)^2 - 2),                                 2*cos(phi/2)^2 + e2^2 - 1,   2*e1*cos(phi/2)*sin(phi/2) - e2*e3*(2*cos(phi/2)^2 - 2)]
% [   2*e2*cos(phi/2)*sin(phi/2) - e1*e3*(2*cos(phi/2)^2 - 2), - 2*e1*cos(phi/2)*sin(phi/2) - e2*e3*(2*cos(phi/2)^2 - 2),                                 2*cos(phi/2)^2 + e3^2 - 1]
 
C = subs(C,e1*sin(phi/2),b1);
C = subs(C,e2*sin(phi/2),b2);
C = subs(C,e3*sin(phi/2),b3);
% substitue (2*cos(phi/2)^2 - 2) expressions in C
% using trig identity cos(phi/2)^2 = 1 - sin(phi/2)^2
%=>
% C =
% [ 2*cos(phi/2)^2 + 2*b1^2 - 1,   2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), - 2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2)]
% [ - 2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e2^2*(2*(1 - sin(phi/2)^2) - 2) - 1,   2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2)]
% [   2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2), - 2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e3^2*(2*(1 - sin(phi/2)^2) - 2) - 1]

C = ...
[ 2*cos(phi/2)^2 + 2*b1^2 - 1,   2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), - 2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2)
 - 2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e2^2*(2*(1 - sin(phi/2)^2) - 2) - 1,   2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2)
 2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2), - 2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e3^2*(2*(1 - sin(phi/2)^2) - 2) - 1];

%Now make substitution from equation 3.97, b0^1+b1^2+b2^2+b3^2 = 1
%=>
% C =
% [ 2*cos(phi/2)^2 + 2*b1^2 - (b0^1+b1^2+b2^2+b3^2),   2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), - 2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2)]
% [ - 2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e2^2*(2*(1 - sin(phi/2)^2) - 2) - (b0^1+b1^2+b2^2+b3^2),   2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2)]
% [   2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2), - 2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e3^2*(2*(1 - sin(phi/2)^2) - 2) - (b0^1+b1^2+b2^2+b3^2)]

C = ...
[ 2*cos(phi/2)^2 + 2*b1^2 - (b0^2+b1^2+b2^2+b3^2),   2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), - 2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2)
 - 2*b3*cos(phi/2) - e1*e2*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e2^2*(2*(1 - sin(phi/2)^2) - 2) - (b0^2+b1^2+b2^2+b3^2),   2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2)
  2*b2*cos(phi/2) - e1*e3*(2*(1 - sin(phi/2)^2) - 2), - 2*b1*cos(phi/2) - e2*e3*(2*(1 - sin(phi/2)^2) - 2), 2*cos(phi/2)^2 - e3^2*(2*(1 - sin(phi/2)^2) - 2) - (b0^2+b1^2+b2^2+b3^2)];

%Now substitue in the given expressions for b in terms of e1,e2, and e3 in
%equation 3.96
C = subs(C,cos(phi/2),b0);
C = subs(C,e1*sin(phi/2),b1);
C = subs(C,e2*sin(phi/2),b2);
C = subs(C,e3*sin(phi/2),b3)
%% 3.23
clear all;close all;clc
%Use eqn 3.98 with eqn 3.102 and equate values
syms b0 b1 b2 b3 b4 bp0 bp1 bp2 bp3 bpp0 bpp1 bpp2 bpp3
%Create the FN matrix in eqn 3.102
FN = [b0^2+b1^2-b2^2-b3^2, 2*(b1*b2+b0*b3), 2*(b1*b3-b0*b2);
        2*(b1*b2-b0*b3), b0^2-b1^2+b2^2-b3^2, 2*(b2*b3+b0*b1);
        2*(b1*b3+b0*b2), 2*(b2*b3-b0*b1), b0^2-b1^2-b2^2+b3^2];
%Create the FB matrix in eqn 3.102    
BN = [bp0^2+bp1^2-bp2^2-bp3^2, 2*(bp1*bp2+bp0*bp3),2*(bp1*bp3-bp0*bp2);
        2*(bp1*bp2-bp0*bp3), bp0^2-bp1^2+bp2^2-bp3^2, 2*(bp2*bp3+bp0*bp1);
        2*(bp1*bp3+bp0*bp2), 2*(bp2*bp3-bp0*bp1), bp0^2-bp1^2-bp2^2+bp3^2];
%Create the BN matrix in eqn 3.102    
FB = [bpp0^2+bpp1^2-bpp2^2-bpp3^2, 2*(bpp1*bpp2+bpp0*bpp3),2*(bpp1*bpp3-bpp0*bpp2);
        2*(bpp1*bpp2-bpp0*bpp3), bpp0^2-bpp1^2+bpp2^2-bpp3^2, 2*(bpp2*bpp3+bpp0*bpp1);
        2*(bpp1*bpp3+bpp0*bpp2), 2*(bpp2*bpp3-bpp0*bpp1), bpp0^2-bpp1^2-bpp2^2+bpp3^2];
    
A = [bpp0 -bpp1 -bpp2 -bpp3;
     bpp1 bpp0 bpp3 -bpp2;
     bpp2 -bpp3 bpp0 bpp1;
     bpp3 bpp2 -bpp1 bpp0];
 
 b = [bp0;bp1;bp2;bp3];
 
eqn3p103 = A*b;
 
 C = FB*BN;
 %Observe how ANS.b0(1) through ANS.b3(1) are equivalent to the four rows
 %of eqn3p103, verifying the composite rotation of equations 3.103 and
 %3.104 using equation 3.98 in 3.102
ANS = solve(FN == C,[b0 b1 b2 b3]);
ANS.b0(1)
ANS.b1(1)
ANS.b2(1)
ANS.b3(1)
eqn3p103