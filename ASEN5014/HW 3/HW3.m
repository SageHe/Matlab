clear all; close all;clc
%% Problem 5.39
% part a -- show that set of vectors given are linearly independent
x1 = [1 2 3]';
x2 = [1 -2 3]';
x3 = [0 1 1]';
%form Gramian out of given vectors
G = [x1'*x1 x1'*x2 x1'*x3;...
        x2'*x1 x2'*x2 x2'*x3;...
        x3'*x1 x3'*x2 x3'*x3];
test = det(G);%since test ~= 0, this set is linearly independent
%part b -- generate orthonormal set using Gram-Schmidt process. First
%results v1 v2 v3 were done by hand, shown in written portion of
%submission. Then tested via Matlab.
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
%% problem 5.40 -- find reciprocal basis set for vector in problem 5.39
%using x1,x2, and x3 from problem 5.39
B = [x1 x2 x3];
%Using info on slide 9 of lecture 11, reciprocal R = ((B'*B)^-1)*B'
R = inv(B'*B)*B';
%From matrix above, reciprocal basis formed from columns of R
r1 = R(:,1);
r2 = R(:,2);
r3 = R(:,3);
%% Problem 5.42 -- Express vector z in terms of basis set xi using reciprocal basis vectors r_i
z = [6 4 -3]';
%From lecture 11 slide 9, using form z = B*a, with RB=I -> R^-1 = B
%-> z = B*a -> R*z = R*B*a = a -> R*z = a
a = R*z;
%%Problem 5.47
y1 = [1 1 1 1]';
y2 = [1 9.9989998e-1 1 1]';
y3 = [-2 -1.9999000 -2 -2]';
%Form matrix A whose columns consist of given y vectors
A = [y1 y2 y3];
%Form Gramian matrix from A
G = A'*A;
%Take determinant of Gramian to determine linear independence
test = det(G);
%det(G) is very nearly 0 -> set of vectors y are nearly linearly dependent,
%and for all intents and purposes are linearly dependent. 
%% Problem 5.49 -- Determine the dimension of the vector space spanned by the given vectors
x1 = [1 2 2 1]';
x2 = [1 0 0 1]';
x3 = [3 4 4 3]';
%By observation, 2*x1 + x2 = x3 -> vectors are not linearly independent ->
%any combo of two of the three vectors form linearly independent set -> 
temp = [x1 x2];
temp1 = [x1 x3];
temp2 = [x2 x3];
G = temp'*temp;
det(G);
G1 = temp1'*temp1;
det(G1);
G2 = temp2'*temp2;
det(G2);
%All determinants of Gramians formed from any two of the three original
%vectors is nonzero -> dimensions of space spanned by set of vectors is
%equal to 2