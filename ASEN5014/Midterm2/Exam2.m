%% Problem 2
% Part a, find basis sets for CS(M),LN(M),RN(M), and RS(M)
M = [1 -3 5 1 6;0 -1 1 0 3;3 -4 10 3 3;1 -1 3 1 0];
y = [1 1 2 -2]';
z = [-1 -2 7 3]';
p = rank(M);
[Q,R] = qr(M);
CS = Q(:,1:p);
LN = Q(:,(p+1):end);
transp_M = M';

[QT,RT] = qr(transp_M);

pt = rank(transp_M);
RS = QT(:,1:pt);
RN = QT(:,(pt+1):end);
% part b, determine components of y in the CS and LN subspaces and
% determine error in a minimum error solution

% y in the CS subspaces using least squares
x_ls_cs = inv(CS'*CS)*CS'*y;
y_ls_cs = CS*x_ls_cs;
% y in the LN subspaces using least squares
x_ls_ln = inv(LN'*LN)*LN'*y;
y_ls_ln = LN*x_ls_ln;
% error in minimum error solution -> error in CS solution
y_err = norm(y - y_ls_cs);
% part c, describe set of all solns to z = Mx for given z, and give minimum
% length soln.

vec = [M z];
vec = rref(vec);

x_z = inv(CS'*CS)*CS'*z;
minlen_z = CS*x_z;

%% Problem 3
% Part a, find all subsets of x in the state space R4 that are subsets
% where solutions x(t) starting 
A = [2 6 0 -6;3 5 0 3;3 -3 4 -1;-3 3 0 5];
[V,D] = eig(A);
sub1 = V(:,1);
sub2 = [V(:,2) V(:,4)];
sub3 = V(:,3);

% part c, find Jordan form for A, along with geometric and algebraic
% multiplicities of each eigenvalue, and determine if A is simple or
% semi-simple 
jordA = jordan(A);