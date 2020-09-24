%% Problem 2, part b
clear all;close all;clc
alpha = [0.6,0.6,1.4];
beta = [1.1,1.75,1];
x = [0 0]';

k = [0:30];

evals = zeros(2,3);
for j = 1:3
    for i = 1:numel(k) - 1
        F = [alpha(j) alpha(j); beta(j)*(alpha(j) - 1) beta(j)*alpha(j)];
        G = [alpha(j);beta(j)*alpha(j)];
        H = [1 1];
        Uk = 1;

        x(:,i+1) = F*x(:,i) + G*Uk;
    end
        figure(j)
        hold on 
        plot(k,x(1,:))
        plot(k,x(2,:))
        evals(:,j) = eig(F)
        title('alpha = 1.4, Beta = 1')
        xlabel('K value');
        ylabel('Magnitude');
%         title('alpha = 0.6, Beta = 1.1 ')
        Y = [H;H*F];
        rank(Y)
end

%% part c
x = [3 1]';

k = [0:30];

evals = zeros(2,3);
for j = 1:3
    for i = 1:numel(k) - 1
        F = [alpha(j) alpha(j); beta(j)*(alpha(j) - 1) beta(j)*alpha(j)];
        G = [alpha(j);beta(j)*alpha(j)];
        H = [1 1];
        Uk = 0;

        x(:,i+1) = F*x(:,i) + G*Uk;
    end
        figure(j+3)
        hold on 
        plot(k,x(1,:))
        plot(k,x(2,:))
        evals(:,j) = eig(F);
        xlabel('K value');
        ylabel('Magnitude');
        title('Alpha = 1.4, Beta = 1, Uk = 0')
        Y = [H;H*F];
        rank(Y)
end
%% part c
x = [3 1]';
y = [];
X0 = [];
for j = 1:3
    O = [];
    for i = 1:10
        F = [alpha(j) alpha(j); beta(j)*(alpha(j) - 1) beta(j)*alpha(j)];
        G = [alpha(j);beta(j)*alpha(j)];
        H = [1 1];
        Uk = 0;

        x(:,i+1) = F*x(:,i) + G*Uk;
        
        O = [O;H*F^(i-1)];
    end
    y = [y;sum(x)];
    x0 = inv(O'*O)*O'*(y(j,[1:end-1]))';
    temp = x0';
    X0 = [X0;temp];
end
Y = [];
for j = 1:3
    O = [];
    F = [alpha(j) alpha(j); beta(j)*(alpha(j) - 1) beta(j)*alpha(j)];
    G = [alpha(j);beta(j)*alpha(j)];
    for i = 1:10
        O = [O;H*F^(i-1)];
    end
    Y(:,j) = O*((X0(j,:))');
end
t = [0:9];
close all
figure
hold on
plot(t,y(1,[1:end-1]))
plot(t,Y(:,1))
title('Predicted vs Real Observations First Alpha Beta Combo')
legend('Real Observations','Predicted Observations')
figure
hold on
plot(t,y(2,[1:end-1]))
plot(t,Y(:,2))
title('Predicted vsReal Observations Second Slpha Beta Combo')
legend('Real Observations','Predicted Observations')
figure
hold on
plot(t,y(3,[1:end-1]))
plot(t,Y(:,3))
title('Predicted vs Real Observation Third Alpha Beta Combo')
legend('Ral Observations','Predicted Observations')
