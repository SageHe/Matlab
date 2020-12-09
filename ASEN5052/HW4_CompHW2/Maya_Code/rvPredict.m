%% Function for Q1
function [rvec,vvec,T,elems] = rvPredict(rvec0,vvec0,t0,mu)
% Inputs - initial pos & vel vectors, initial time, mu, t?
% Outputs - r and v vectors at time t
r0 = norm(rvec0);

%% Step 1 - Orb Elem from r0Vec v0Vec

% hVec
hvec = cross(rvec0,vvec0);
h = norm(hvec);
hhat = hvec/h;

% e
evec = (1/mu)*cross(vvec0,hvec) - rvec0/r0;
e = norm(evec);
ehat = evec/e;

% a,T
a = (h^2)/mu/(1-(e^2));
T = sqrt(a^3/mu)*2*pi;

% i 
i = acos(hhat(3)); %rad
if i < 0  %making sure 0<i<pi
    i = i+pi;
elseif i > pi
    i = i-pi;
end

% Omega
nOmega = cross([0 0 1],hhat);
nOmegahat = nOmega/norm(nOmega);
Omega = atan(nOmegahat(2)/nOmegahat(1)); %rad
%make sure 0<=Omega<2pi

% womega
nOmegahatPerp = [-cos(i)*sin(Omega) cos(i)*cos(Omega) sin(i)];
womega = atan(dot(ehat,nOmegahatPerp)/dot(ehat,nOmegahat));

%% Step 2 - Solve Keplers for t, to find f
% solve for f as an initial guess for E...
ehatPerp = [-sin(womega)*cos(Omega) -sin(womega)*sin(Omega) 0];
Eg = atan2(dot(rvec0,ehatPerp),dot(rvec0,ehat));
E = newtons(e,t0,mu,a,Eg);
% now solve for f using E
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

%% Step 3 - convert alpha to rvec, vvec

% new rvec:
rvec = (a*(1-(e^2))/(1+e*cos(f)))*((cos(f)*ehat) + sin(f)*cross(hhat,ehat));
r = norm(rvec);
% new vvec using Eqs. 3.52,3.53,3.58,3.68,3.69
if e==1 %parabola
    X = h*(tan(f/2)-tan(Eg/2))/sqrt(mu); % assumes t0 @ peri and f0 = 0
else % ellipse or circle, assuming no hyperbolas
    X = sqrt(a)*(E-Eg); % assumes t0 @ peri and E0 = 0
end

alpha = 1/a;
z = alpha*(X^2);

if z>0
    S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
    C = (1-cos(sqrt(z)))/z;
elseif z<0
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    C = (cosh(sqrt(-z))-1)/-z;
else
    S = 1/6;
    C = 1/2;
end

f_dot = sqrt(mu)*(alpha*(X^3)*S - X)/r/r0;
g_dot = 1 - ((X^2)*C/r);

vvec = f_dot*rvec0 + g_dot*vvec0;

v0 = norm(vvec0);
Energy = (.5*v0^2) + (.5*mu^2/(r0^2)) - (mu/r0);
elems = [a e i womega Omega E];




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



