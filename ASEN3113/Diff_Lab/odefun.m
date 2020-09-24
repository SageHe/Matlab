function [yPrime] = odefun(t,y)
%t = input time
%y = input vector, y = [y1,y2]
%yPrime = vector of outputs yPrime = [y1Prime, y2Prime]

yPrime = zeros(2,1);
yPrime(1) = y(2);
yPrime(2) = -y(1);
end