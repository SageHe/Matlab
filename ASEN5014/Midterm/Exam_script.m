%% Problem 1
% Part a
M = [1 0 0;0 1 0;1 0 1];
det(M'*M)
rank(M)
% Part b
clear M
M = [1 0 -1;1 2 1];
rank(M)
%det(M'*M)
% Part j
M1 = [1 0 1;0 0 1;0 1 1;1 0 0];
M2 = [1 0 0;0 1 0;1 0 1];
rank(M1)
rank(M2)
%% Problem 2
% Part a
B = [3 1 2;1 1 0;-1 2 3;2 -1 -2;4 -2 4;2 3 1]; 
C = [5 2 0;1 0 2;2 -3 1;0 3 0;8 6 -8;3 -1 5];
z = [18 2.6 3.1 3.0 34 7.1]';
RB = inv(B'*B)*B';
RC = inv(C'*C)*C';
zb = RB*z;
zc = RC*z;
% Part b
w = [3 24 3 12 1 6];
test = [B w'];
det(test'*test)
rank(test)
% Part c -- test all vectors in B and C to determine if every vector in C
% is linearly independent of those in B and vice versa
rank([B C(:,1)])
rank([B C(:,2)])
rank([B C(:,3)])
rank([C B(:,1)])
rank([C B(:,2)])
rank([C B(:,3)])
% Part d
T = RB*C