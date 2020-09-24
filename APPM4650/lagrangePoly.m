%Ryan Cameron
%APPM 4650
%Homework 4
%This function finds the lagrange polynomial for the given x-points and
%evaluates it at the desired number

function [poly,eval] = lagrangePoly(xPoints,yPoints,xReal)
syms x
deg = length(xPoints)-1;
poly = 0;
for i = 1:length(xPoints) %Number of lagrange polynomials that make up the whole
    %Create the numerator and denominator, then divide
    num = 1;
    den = 1;
    for j = 1:length(xPoints)
        if j~=i
            new = (x-xPoints(j));
            num = num*new;
            
            new = (xPoints(i) - xPoints(j));
            den = den * new;
        end
    end
    newL = num/den;
    poly = poly + (newL*yPoints(i));
end
eval = subs(poly,x,xReal);
eval = double(eval);
poly = vpa(poly,5);
end