clear all; close all;
%This script calculates the binomial coefficient without using the built-in
%matlab factorial function, with n greater than k, n and k are nonnegative
%integers

Nfact = 1;
Kfact = 1;
NKfact = 1;
N = input('Input the desired N value of the binomial coefficient\n');
K = input('Input the desired K value of the binomial coefficient\n');
 
 while N < K
     fprintf('DANGER! N must be greater than K \n');
     K = input('Please input a value of K that is less than the value of K \n');
 end
 
 while (N < 0 || K < 0)
     fprintf('DANGER! N and K must be positive integers \n');
     while (N < 0)
         N = input('Please enter a positive value of N \n');
     end
     while (K < 0)
         K = input('Please enter a positive value of K \n');
     end
 end

 for i = 1:N
     Nfact = Nfact * i;
 end
 
 for i = 1:K
     Kfact = Kfact * i;
 end
     
 for i = 1:(N - K)
     NKfact = NKfact * i;
 end
 
 N_choose_K = Nfact/((NKfact).*Kfact);
 
 fprintf('%d choose %d is equal to %d \n',N,K,N_choose_K);