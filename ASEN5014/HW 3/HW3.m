clear all; close all;clc
%% Problem 5.39
% part a - show that set of vectors given are linearly independent
x1 = [1 2 3]';
x2 = [1 -2 3]';
x3 = [0 1 1]';
G = [x1'*x1 x1'*x2 x1'*x3;...
        x2'*x1 x2'*x2 x2'*x3;...
        x3'*x1 x3'*x2 x3'*x3];
test = det(G);%since test ~= 0, this set is linearly independent
%part b - generate orthonormal set using Gram-Schmidt process
v1 = [1 2 3]';
v2 = [4/7 -20/7 12/7]';
v3 = [([0 1 1] - [5/14 10/14 15/14] - [-1/10 -2/10 -3/10])'];
vhat1 = v1/norm(v1);
vhat2 = v2/norm(v2);
vhat3 = v3/norm(v3);
%testing with ops. done in code
v1 = x1;
v2 = x2 - (dot(v1,x2)/dot(v1,v1))*v1;
v3 = x3 - (dot(v1,x3)/dot(v1,v1))*v1 - (dot(v2,x3)/dot(v2,v2))*v2;
vhat1 = v1/norm(v1);
vhat2 = v2/norm(v2);
vhat3 = v3/norm(v3);