function dydt = revenge_of_the_szr(t, y, alpha, beta, gamma, N0)
%Code below was given in project3 pdf
% Evaluates the right hand side of the SZR model for equations S(t) and Z(t).
% Here, y(t) = [S(t); Z(t)], so y(1) = S(t) and y(2) = Z(t).
dSdt = -beta*y(1)*y(2);
if t < 14
    dZdt = beta*y(1)*y(2) + gamma*(N0 - y(1) - y(2)) - alpha*y(1)*y(2);
else
    dZdt = beta*y(1)*y(2) + gamma*(N0 - y(1) - y(2)) - alpha*y(1)*y(2) - y(1)*y(2)*.00005*t;
end

dydt=[dSdt;dZdt];
end