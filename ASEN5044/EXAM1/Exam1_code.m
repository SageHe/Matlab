clear; close all;clc
%Define given constants in order to convert CT LTI to DT LTI
g = 9.81; %m/s^2
l = 1; %meters
m = 1; % Kg
M = 2; %Kg
deltat = 0.05; %seconds
%Define A and B matries to form augmented matrix in order to determin F and
%G matrices using matrix exponential (expm function)
A = [0 1 0 0;0 0 (m*g)/M 0;0 0 0 1;0 0 ((g/l)/(1 - (m/(M+m)))) 0];

B = [0; 1/M; 0; 1/(M*l)];

C = [1 0 -l 0];
D = 0;
Ahat = [A B];
Ahat = [Ahat;zeros(1,5)];

%Use matrix exponential to compute F and G
matexp = expm(Ahat*deltat);
%Pull F and G matrices out of matexp and define
F = matexp([1:4],[1:4])
G = matexp([1:4],end)
H = C
M = D
%Determine stability of system by observing eigenvalues of F
[V,D] = eig(F);
D
%Since there is an eigenvalue greater than 1, the system is asymptotically
%unstable 

%Form O matrix to determin observability of the system
O = [H;H*F;F*F^2;H*F^3];

rank(O)

%Since the rank is equal to the number of measurements, which is 4, the
%system is ovservable 


