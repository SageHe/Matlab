%%Question 1
%Part a
T1 = [-10000 -10000 10000];
T2 = [-10000 10000 10000];
T3 = [10000 -10000 10000];
T4 = [10000 10000 10000];
rangeT1 = sqrt(T1(1)^2 + T1(2)^2 + T1(3)^2);
azT1 = atan2d(T1(1),T1(2));
elT1 = asind(T1(3)/rangeT1);

rangeT2 = sqrt(T2(1)^2 + T2(2)^2 + T2(3)^2);
azT2 = atan2d(T2(1),T2(2));
elT2 = asind(T2(3)/rangeT2);

rangeT3 = sqrt(T3(1)^2 + T3(2)^2 + T3(3)^2);
azT3 = atan2d(T3(1),T3(2));
elT3 = asind(T3(3)/rangeT3);

rangeT4 = sqrt(T4(1)^2 + T4(2)^2 + T4(3)^2);
azT4 = atan2d(T4(1),T4(2));
elT4 = asind(T4(3)/rangeT4);
%Part b
G = [10000/17320.6 10000/17320.6 1;...
    10000/17320.6 -10000/17320.6 1;...
    -10000/17320.6 +10000/17320.6 1;...
    -10000/17320.6 -10000/17320.6 1];
dp = [33.3 32.2 32.1 33.5]';
%Part c
H = inv(G'*G);
dops = sqrt(diag(H));
Edop = dops(1);
Ndop = dops(2);
Tdop = dops(3);
%Part d
lssoln = inv(G'*G)*G'*dp;
%Recalc if observer is off by 1000m in y
G = [10000/17320.6 11000/17320.6 1;...
    10000/17320.6 -9000/17320.6 1;...
    -10000/17320.6 +11000/17320.6 1;...
    -10000/17320.6 -9000/17320.6 1];
dp = [17353.9-16763.5 17352.8-17916.4 17352.7-16763.5 17354.1-17916.4]';
lssoln_shift = inv(G'*G)*G'*dp;
lssoln_shift = [0 1000 0]' + lssoln_shift;