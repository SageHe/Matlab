function [astar] = lamb(r0,rf,TOF,mu)
dtheta_star = acos((dot(r0,rf))/(norm(r0)*norm(rf)));

dtheta_star = 2*pi - dtheta_star;

c = sqrt(norm(r0)^2 + norm(rf)^2 - 2*norm(r0)*norm(rf)*cos(dtheta_star));
s = .5*(norm(r0) + norm(rf) + c);

TOF_p = (1/3)*sqrt(2/mu)*(s^(3/2) + (s - c)^(3/2));

am = (s/2);
nm = sqrt(mu/am^3);
alpha_m = pi;
beta_m0 = 2*asin(sqrt((s - c)/s));
beta_m = -beta_m0;

a = am;

alpha0 = @(a) 2*asin(sqrt(s/(2*a)));
beta0 = @(a) 2*asin(sqrt((s-c)/(2*a)));

alpha = @(a) 2*pi - 2*asin(sqrt(s/(2*a)));
beta = @(a) -(2*asin(sqrt((s-c)/(2*a))));

n = @(a) sqrt(mu/a^3);

fun = @(a) TOF - (1/n(a))*(alpha(a) - beta(a) - (sin(alpha(a)) - sin(beta(a))));

astar = fsolve(fun,am);
end