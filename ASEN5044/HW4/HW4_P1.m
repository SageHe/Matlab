clear all;close all;clc

A = [0 1 0 0;-2 0 1 0;0 0 0 1;1 0 -2 0];
B = [0 0; -1 0; 0 0; 1 1];

A_aug = [A B];
pad = zeros(2,6);
A_aug = [A_aug;pad];

A_hat = expm(A_aug*.05);

F = A_hat([1:4],[1:4]);
G = A_hat([1:4],[5:6]);

[V,D] = eig(A);
%% Part b 
H = [1 0 0 0;0 1 0 -1];
Y = [H;H*F;H*F^2;H*F^3];

rank(Y)
%% Part c/d
load('hw4problem1data.mat');
Udata(end,:) = [];
O = [];
% theta = 0;
for i = 1:size(Ydata,1)
    theta = 0;
    temp = H*F^i;
    O = [O;temp];
    for j = 1:i
        theta = theta + F^(i - j)*(G*Udata(j,:)');
    end
    theta = H*theta;
    phi(i,:) = Ydata(i,:) - theta';
end
phi = reshape(phi',[200,1]);
x0 = inv(O'*O)*O'*phi;

for i = 1:100
    y(:,i) = H*F^(i-1)*x0;
end
t = linspace(0,5,size(Ydata,1));
figure
hold on
plot(t,y(1,:))
plot(t,y(2,:))
title('Predicted Positions')
xlabel('Time (s)')
ylabel('Position (m)')
hold off
legend('Predicted Position 1','Predicted Position 2')
figure
hold on
plot(t,Ydata(:,1))
plot(t,Ydata(:,2))
title('Recorded Measurements')
xlabel('Time (s)')
ylabel('Position (m)')
hold off
legend('Calculated Position 1','Calculated Position 2')
figure
hold on
plot(t,abs(y(1,:)-Ydata(:,1)'))
plot(t,abs(y(2,:)-Ydata(:,2)'))
title('Error Between Predicted and Recorded Measurements')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Difference in Predicted vs Calculated Position 1','Difference in Predicted vs Calculated Position 2')