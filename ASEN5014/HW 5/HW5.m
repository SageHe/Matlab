%% Problem 6.36 from textbook
% Find all nontrivial solutions of Ax = 0, i.e. null space of A
A = [26 17 8 39 35;17 13 9 29 28;8 9 10 19 21;39 29 19 65 62;35 28 21 62 61];
[Q,R] = qr(A);
%% Problem 6.41 from textbook
% Find least squares soluttion to a and b
y = [5 1 7]';
M = [1 2;1 -2;1 5];
temp = inv(M'*M)*M';
abvec = temp*y;
%% Problem 6.45 from textbook
%Solve for C and alpha and predict time for one mile
T = [185 79.6 60 37.9 11.5]';
D = [26.2 12.4 9.5 6.2 2]';

y = log(T);
M = [ones(5,1) log(D)];

calphavec = inv(M'*M)*M'*y;