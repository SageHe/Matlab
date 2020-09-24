clear all;close all

h = input('What is the desired step size you wish to use? \n');
% t0 = input('Input the starting time \n');
% t1 = input('Input the ending time \n');

%t = t0:h:t1;
y = [];
y(1) = input('What is the initial condition (Y(0) = ? )? \n');
p = input('What is the monthly payment amount? \n');
r = input('What is the interest rate on the loan? \n');
n = 1;
m = 0;
% for n = 1:(length(t) - 1)
%     y = [y (y(n) + h*yprime(t(n),y(n)))];
% end
while y > 0
    y = [y (y(n) + h*yprime1(y(n),r,p))];
    n = n + 1;
    m = m + h;
end
t= 0:h:(n - 1)*h;
y_exact = ((exp(r*t))*(y(1)*r - 12*p) + 12*p)/r;   
hold on

plot(t,y)
plot(t,y_exact)
legend('Numerical Solution','True Solution')
title('Analytical vs. True solution')
xlabel('t(years)')
ylabel('Y(dollars)')

