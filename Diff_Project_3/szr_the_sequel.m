function dydt = szr_the_sequel(t, y, alpha, beta, gamma, N0, rho)
%Code below was given in project3 pdf
% Evaluates the right hand side of the SZR model for equations S(t) and Z(t).
% Here, y(t) = [S(t); Z(t)], so y(1) = S(t) and y(2) = Z(t).
dSdt = -beta*y(1)*y(2) + rho*y(2);
dZdt = beta*y(1)*y(2) + gamma*(N0 - y(1) - y(2)) - alpha*y(1)*y(2) - rho*y(2);
dydt=[dSdt;dZdt];
end