gamma = 1.4;

eqn1 = @(M) ((1/M)*((2/(gamma + 1))*(1 + ((gamma - 1)/2)*M^2))^((gamma + 1)/(2*(gamma - 1)))*((((gamma + 1)/2)*M^2)/(1 + ((gamma - 1)/2)*M^2))^(gamma/(gamma - 1))*(((2*gamma)/(gamma + 1))*M^2 - ((gamma - 1)/(gamma + 1)))^(-1/(gamma - 1))) - 1.2;

soln = fzero(eqn1,[1,2.5]);