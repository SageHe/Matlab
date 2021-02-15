Ephem = Meeus(2456300);

a = Ephem.Mars.a;
e = Ephem.Mars.e;
i = Ephem.Mars.i;
w = wrapTo360(Ephem.Mars.w);
Omega = wrapTo360(Ephem.Mars.Omega);
nu = wrapTo360(Ephem.Mars.nu);
M = wrapTo360(Ephem.Mars.M);
% mu = 3.986004415e5;
mu = 1.32712440018e11;
X = kep2car([a e i Omega w nu M],mu,'degrees');
[R,V] = calcposvel(a,e,i,Omega,w,nu,mu);
% [R2,V2] = kep2cart([a e deg2rad(i) Omega deg2rad(w) nu M],mu);