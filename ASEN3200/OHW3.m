%% Question 3, 4.18
Omega = deg2rad(75);
i = deg2rad(50);
w = deg2rad(80);

Q = [-sin(Omega)*cos(i)*sin(w)+cos(Omega)*cos(w) -sin(Omega)*cos(i)*cos(w)-cos(Omega)*sin(w) sin(Omega)*sin(i);
    cos(Omega)*cos(i)*sin(w)+sin(Omega)*cos(w) cos(Omega)*cos(i)*cos(w)-sin(Omega)*sin(w) -cos(Omega)*sin(i);
    sin(i)*sin(w) sin(i)*cos(w) cos(i)];

r = [6578 0 0]';

Q*r

v = [0 11.546 0]';

Q*v
%% Question 4
r = [4.973e3 -1.798e3 -1.748e3];
h = [-3110.39 -9623.7 1048.54];
v = [-.5713 .4174 2.136];

e = (1/22030)*cross(v,h) - r/norm(r)

