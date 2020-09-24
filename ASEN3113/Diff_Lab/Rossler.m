function [yPrime] = Rossler(t,y)
%t = input time
%y = input vector, y = [y1,y2,y3]
%yPrime = vector of outputs yPrime = [y1Prime, y2Prime,y3prime]

a = 0.1;
b = 0.1;
c = 14;
yPrime = zeros(3,1);
yPrime(1) = -y(2) - y(3);
yPrime(2) = y(1) + a*y(2);
yPrime(3) = b + y(3)*(y(1) - c);
end