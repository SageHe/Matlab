function [yPrime] = lorentz(t,y)
%t = input time
%y = input vector, y = [y1,y2,y3]
%yPrime = vector of outputs yPrime = [y1Prime, y2Prime,y3prime]

sigma = 10;
rho = 28;
beta = 8/3;
yPrime = zeros(3,1);
yPrime(1) =sigma*(y(2) - y(1));
yPrime(2) = y(1)*(rho - y(3)) - y(2);
yPrime(3) = y(1)*y(2) - beta*y(3);
end