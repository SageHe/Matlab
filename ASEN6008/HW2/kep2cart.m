function [rvec,vvec] = kep2cart(orb_el,mu)
a = orb_el(1);
e = orb_el(2);
inc = orb_el(3);
RAAN = deg2rad(orb_el(4));
w = orb_el(5);
ths = (orb_el(6));
M = orb_el(7);
E = mean2eccentric(M,e);
r = a*(1-e^2)/(1+e*cos(ths));
V = sqrt(2*mu/r - mu/a);
fpa = e*sin(ths)/(1+e*cos(ths));
rrotvec = [r, 0, 0];
vr = V*sin(fpa); % km/s
vth = V*cos(fpa); % km/s
vrotvec = [vr, vth, 0];
th = ths + w;
C = [cos(RAAN)*cos(th)-sin(RAAN)*cos(inc)*sin(th),...
    -cos(RAAN)*sin(th)-sin(RAAN)*cos(inc)*cos(th), sin(RAAN)*sin(inc);...
    sin(RAAN)*cos(th)+cos(RAAN)*cos(inc)*sin(th),...
    -sin(RAAN)*sin(th)+cos(RAAN)*cos(inc)*cos(th), -cos(RAAN)*sin(inc);...
    sin(inc)*sin(th), sin(inc)*cos(th), cos(inc)];
rvec = C*rrotvec';
vvec = C*vrotvec';
end

