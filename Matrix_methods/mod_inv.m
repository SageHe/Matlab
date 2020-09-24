%clear all; close all; clc

a = input('what number do you want to find the inverse of?\n');
b = input('of what modulo?\n');
%Check if modular inverse exists
if gcd(a,b) ~=1
    error('There is no modular inverse of this integer for the specified modulo')
end

while a < 0
    a = a + b;
end

c = 1;
c = c + b;

while mod(c,a) ~= 0
    c = c + b;
end

modular_inverse = c/a;

fprintf('the inverse of %d modulo %d is %d\n',a,b,modular_inverse)