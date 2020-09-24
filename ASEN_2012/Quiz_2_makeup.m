clear all;close all;clc
%Sage Herrin, ASEN 2012 2PM lecture
%Data created: 11/3/17
%Due: 11/3/17
%Purpose: Determine weighting matrix for situation, compute coefficients of
%line of best fit, determien values of k and m and their uncertainties.
%Given:T and M values and uncertainty in T
%Required:Q matrix and partia derivatives based on equation 3.48 of John R.
%Taylor's Intro to Error Analysis
%Assumptions: Error in T is random.

M = [55:50:455]; %Given values of M from the problem
T = [.496 .645 .761 .867 .957 1.037 1.113 1.194 1.254];

%To build weight matrix, use given uncertinty in T value of 0.0005
T_uncert = 0.0005;
W = eye(9);
W = W*(1/(((T_uncert)^2)^2)); %Weight matrix

%Compute coefficients of line of best fit using matrix multiplication
M = M(:);
T = T(:);
T_sqr = T.^2;
M = [M ones(9,1)];

P = M\T_sqr;

alpha = P(1); %slope coefficient 
Beta = P(2); %intercept coefficient 

Y = alpha*M + Beta;

hold on
plot(M,T_sqr)
plot(M,Y)
title('Mass Spring System Oscillation')
xlabel('Mass (g)')
ylabel('Period T (sec)')
legend('Mass vs Time^2','Weighted least squares best fit line')
%Determine values of m and k using given equation T^2 = ((4pi^2)/k)*M +
%((4pi^2*m)/k) = alphaM + beta, implying alpha = (4pi^2)/k and beta =
%(4(pi^2)m)/k
%Start by solving for k. Manipulating equation gives us 4pi^2/alpha = k
k = (4*(pi^2))/alpha;
%Now k can be used to calculated m by manipulating given equation
m = (Beta * k)/(4 * (pi^2));

%Calculate uncertainties in k and m
%This can be done by using the general form for uncertainties using partial
%derivatives using alpha = (4pi^2)/k and beta = (4(pi^2)m)/k and
%uncertainties in alpha and beta using Q matrix = inv(transpose(A)*W*A)

Q = inv(transpose(M)*W*M);

alpha_part_WRT_k = abs((-4*(pi^2))/k^2);
Beta_part_WRT_k = abs((-4*(pi^2)*m)/k^2);
Beta_part_WRT_m = abs(4*(pi^2)/k);
%Use principle of uncertainty in alpha equal to partial of alpha with
%respect to k multiplied by the uncertainty in k to solve for uncertainty
%in k, use same principle to solve for uncertainty in m as well.
k_uncert = Q(1,1)/alpha_part_WRT_k;
m_uncert = ((Q(2,2) - (Beta_part_WRT_k*k_uncert)) / (Beta_part_WRT_m));



