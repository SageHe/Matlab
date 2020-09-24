mdot = 100;
A = 0.56;
T_t = 288.8;
P_t = 101325;
g = 1.4;
R = 287;

eqn = @(M) sqrt(g/R)*M*(1 + ((g - 1)/2)*M^2)^(-((g + 1)/(2*(g - 1)))) - (mdot/A)*(sqrt(T_t)/P_t);

soln = fzero(eqn,[1,2]);