clear all;close all
%Create the A matrix and b matrix and solve using both Prof. Jackson's
%metod and text book equations
A = [1;2;3;4;5;6];
one = ones(length(A),1);
d = [5;14.4;23.1;32.3;41.0;50.4];
A = [A one];
x = 1:6;
P = inv((transpose(A)*A))*transpose(A)*d;

b = P(2);

m = P(1);

lambda = m * 2;
%Taylor book equations used to calculate A,B,and delta
y = b + m*x;
delta = 6*sum(x.^2) - (sum(x)).^2;
newa = ((sum(x.^2)*(sum(y)) - sum(x)*sum(x.*y))/delta);
newb = ((6*sum(x.*y) - sum(x)*sum(y))/delta);
%Calculate uncertainty in y
sigmay = [];
for i = 1: length(d)
    sigmay = [sigmay (d(i) - newa - newb*x(i)).^2];
end
sigmay = sum(sigmay);
sigmay = sqrt((1/4)*sigmay);
sigmab = sigmay * sqrt(6/delta);
sigmab = sigmab*2;
%plot points given and best straight-fit line
hold on
plot(x,d,'*')
plot(x,y)

title('Kundts Tube')
xlabel('Node Number')
ylabel('Position (cm)')
legend('plotted data','straight-line fit')

fileID = fopen('leastsquares8.15.txt','w');
fprintf(fileID,'The wavelenth was calculated as %.5f with an uncertainty of %.2f \n',lambda,sigmab);
fprintf(fileID,'The straight-line equation was calculated as 9.0286X - 3.9');
fclose(fileID);
