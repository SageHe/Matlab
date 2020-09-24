function dydt = szr(t, y, alpha, beta, gamma, N0)
%Code below was given in project3 pdf
% Evaluates the right hand side of the SZR model for equations S(t) and Z(t).
% Here, y(t) = [S(t); Z(t)], so y(1) = S(t) and y(2) = Z(t).
dydt = [-beta*y(1)*y(2);
beta*y(1)*y(2) + gamma*(N0 - y(1) - y(2)) - alpha*y(1)*y(2)];
end