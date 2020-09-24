clear;close all;clc
%% Problem 1, part b converting CT to DT stochastic LTI
A = [0 1 ;-100 -10];
B = zeros(2);
C = [1 0;0 0.2];
D = 0;

dt = 0.2;

sys = ss(A,B,C,D);

sys = c2d(sys,dt);

[F,G,H,M] = ssdata(sys);
% Compute Q matrix
gamma = [0 1]';
W = 10;

Z = [-A gamma*W*gamma';zeros(2) A']
Z = dt*Z;

ez = expm(Z);

Q = ez([3:4],[3:4])'*ez([1:2],[3:4]);
%% Problem 4
%part a
clear;close all;clc
%Define given covariance matrix R
R = [8 5.15 6.5;5.15 5 -4.07;6.5 -4.07 50];
%Find Sv using cholesky factorization of R
Sv = chol(R,'lower');
%define qk for random samples and compute yk for 100 measurements using 
%yk = mk +Sv*qk, where mk is desired mean (1 in this case)
for i = 1:100
    qk(:,i) = randn(1,3)';
    yk(:,i) = 1 + Sv*qk(:,i);
end
figure(1)
scatter(yk(1,:),yk(2,:))
xlim([-6 10])
ylim([-20 25])
xlabel('y1 (m)')
ylabel('y2 (m)')
title('y1 vs y2')
figure(2)
scatter(yk(1,:),yk(3,:))
xlim([-6 10])
ylim([-20 25])
xlabel('y1 (m)')
ylabel('y3 (m)')
title('y1 vs y3')
figure(3)
scatter(yk(2,:),yk(3,:))
xlim([-6 10])
ylim([-20 25])
xlabel('y2 (m)')
ylabel('y3 (m)')
title('y2 vs y3')
%part b
%Using the built-in covariance function "cov", the sample covariance matrix
%can be determined 
covmat = cov(yk')
%part c
%For use of first 3 observations
ind = 1;
for T = [3 10 100]
    % H = ones(3*T,3);
    H = [];
    for i = 1:T
        H = [H;eye(3)];
    end
    for i = 1:3:3*T
        bigR(i:i+2,i:i+2) = R;
    end
    yk = reshape(yk',300,1);
    xLS(:,ind) = inv(H'*inv(bigR)*H)*H'*inv(bigR)*yk(1:3*T);
    %determine error cov. matrix 
    qk = reshape(qk',300,1);
    els(:,ind) = -inv(H'*inv(bigR)*H)*(H'*inv(bigR))*qk(1:3*T);
    ind = ind+1;
end
xLS(:,1) % x estimate using first 3 measurements
xLS(:,2) % x estimate using first 10 measurements
xLS(:,3) % x estimate using 100 measurements
%respective estimation error cov. matrices
els(:,1)
els(:,2)
els(:,3)
%part d 
data = load('hw6problem5data.csv');
%construct H and R to have proper dimensions 
clear bigR
H = [];
for i = 1:30
    H = [H;eye(3)];
end
for i = 1:3:90
    bigR(i:i+2,i:i+2) = R;
end
data = reshape(data,90,1);
xLS = inv(H'*inv(bigR)*H)*H'*inv(bigR)*data;
els = -inv(H'*inv(bigR)*H)*(H'*inv(bigR))*qk(1:90);
%part e
%using unweighted least squares is the same process as in part d but
%without the use of the R matrix
xLS = inv(H'*H)*H'*data;
els = -inv(H'*H)*(H'*qk(1:90));
%part f
%RLLS 

