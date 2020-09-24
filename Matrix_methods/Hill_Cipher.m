clear all; close all; clc
% Mod. inverses - check multiplicative inverse of a modulo m exists
% => gcd(a,m) = 1. Solve equation ar + ms = 1, r is modular inverse of 
% a modulo m

%Creates message string
message = input('What is the message you would like to encrypt\n');
%Turns string into numbers
message = double(message);
%Reduces string by 97 to reduce determinant magnitude
message = message - 97;
% Random nxn encryption mat., n=string length
gen_key;






