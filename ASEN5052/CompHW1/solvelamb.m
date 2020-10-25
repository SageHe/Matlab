function [astar,V,tp] = solvelamb(r1,r2,tspan,mu)
c = norm(r2 - r1);
s = .5*(norm(r1) + norm(r2) + c);
theta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
tp = (1/sqrt(mu))*(sqrt(2)/3)*((s^(3/2)) - sign(sin(theta))*(s - c)^(3/2));
beta_m = 2*asin(sqrt((s - c)/s));
if (theta >= pi) && (theta <= 2*pi)
    beta_m = -beta_m;
end
tm = (1/sqrt(mu))*sqrt(s^3/8)*(pi - beta_m + sin(beta_m));
syms a 
alpha = 2*asin(sqrt(s/(2*a)));
beta = 2*asin(sqrt((s-c)/(2*a)));
if (theta >= pi)
    beta = -beta;
end
if (tspan(2) - tspan(1) > tm)
    alpha = 2*pi - alpha;
end
f = sqrt(mu)*(tspan(2) - tspan(1)) - a^(3/2)*(alpha - beta - (sin(alpha) - sin(beta)));
fprime = diff(f);
tol = 0.0001;
astar = (s/2)+1;
while abs(double(subs(f,astar))) > tol
    astar = astar - double(subs(f,astar))/double(subs(fprime,astar));
end
A = sqrt(mu/(4*astar))*cot(double(subs(alpha,astar))/2);
B = sqrt(mu/(4*astar))*cot(double(subs(beta,astar))/2);
u1 = r1/norm(r1);
u2 = r2/norm(r2);
uc = (r2 - r1)/norm(r2 - r1);
V = (B + A)*uc + (B - A)*u1;
%part i -- min energy transfer ellipse converged to a semi-major axis of
%1.309 with a transfer time of 4.5885 seconds or 
beta_m = -beta_m;
tmflip = (1/sqrt(mu))*sqrt(s^3/8)*(pi - beta_m + sin(beta_m));
%part ii
min_e = (norm(r2) - norm(r1))/c;
min_e_a = (norm(r1) + norm(r2))/2;
alpha = 2*asin(sqrt(s/(2*min_e_a)));
beta = 2*asin(sqrt((s-c)/(2*min_e_a)));
tmin_e  = (min_e_a^(3/2)*(alpha - beta - (sin(alpha) - sin(beta))))/sqrt(mu);
beta = -beta;
alpha = 2*pi - alpha;
tmin_e_flip = (min_e_a^(3/2)*(alpha - beta - (sin(alpha) - sin(beta))))/sqrt(mu);
%part iii -- V = [0.4829 1.0115 0]

end