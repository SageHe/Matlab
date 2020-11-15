%% ASEN 5519 
% function for 1ii, elements to r,v
function [elemRvec, elemVvec] = rvElements(a,e,i,womega,Omega,t,mu)

E = newtons(e,t,mu,a,1);
% now solve for f using E
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

h = sqrt(a*mu*(1-e^2));

rxBar = (h^2)/mu/(1+e*cos(f))*[cos(f) sin(f) 0];
vxBar = mu/h*[-sin(f) e+cos(f) 0];

QxX = [-sin(Omega)*cos(i)*sin(womega)+cos(Omega)*cos(womega), -sin(Omega)*cos(i)*cos(womega)-...
    cos(Omega)*sin(womega), sin(Omega)*sin(i); cos(Omega)*cos(i)*sin(womega)+sin(Omega)*cos(womega),...
    cos(Omega)*cos(i)*cos(womega)-sin(Omega)*sin(womega), -cos(Omega)*sin(i); sin(i)*sin(womega), ...
    sin(i)*cos(womega), cos(i)];

elemRvec = QxX*rxBar';
elemVvec = QxX*vxBar';


%% Newton's Method

function E = newtons(e,t,mu,a,Eg)

% LHS
n = sqrt(mu/(a^3));
M = n*t;

% inital conditions
fg = Eg - e*sin(Eg);
check = fg - M;
E = Eg;

while abs(check) > 0.00001
    mstr = E - e*sin(E) - M;
    mstr_prime = 1 - e*cos(E);
    E = E - (mstr/mstr_prime);
    check = E - e*sin(E) - M;
end

end

end
