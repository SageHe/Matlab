%% Problem 5.41 -- express vector z in terms of orthonormal basis set vi of problem 5.39
%recreating the orthonormal basis set
x1 = [1 2 3]';
x2 = [1 -2 3]';
x3 = [0 1 1]';
v1 = x1;
v2 = x2 - (dot(v1,x2)/dot(v1,v1))*v1;
v3 = x3 - (dot(v1,x3)/dot(v1,v1))*v1 - (dot(v2,x3)/dot(v2,v2))*v2;
vhat1 = v1/norm(v1);
vhat2 = v2/norm(v2);
vhat3 = v3/norm(v3);
%use vhats to express Z in terms of orthonormal basis set
z = [6 4 -3]';
V = [vhat1 vhat2 vhat3];
x = inv(V)*z;
%% Problem 3
V = [2 1;0 -1;1 1;2 2;2 1];
[Q,R] = qr(V);
orth_comp = Q(:,3:5);