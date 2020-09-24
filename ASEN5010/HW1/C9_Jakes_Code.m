N = 10000;
t = linspace(0,42,N);
tstep = t(2)-t(1);
psi = 40*pi/180;
theta = 30*pi/180;
phi = 80*pi/180;

for i = 1:N
omega = 20*pi/180*[sin(0.1*t(i)) 0.01 cos(0.1*t(i))]';
[dpsi,dtheta,dphi] = kde(psi(i),theta(i),phi(i),omega);

psi(i+1) = psi(i)+dpsi*tstep;
theta(i+1) = theta(i)+dtheta*tstep;
phi(i+1) = phi(i) + dphi*tstep;
end

disp(norm([phi(end),theta(end),psi(end)]))
plot(psi)
hold on
plot(theta)
plot(phi)


function [dpsi,dtheta,dphi] = kde(psi,theta,phi,omega)
c1=cos(psi);
c2=cos(theta);
c3=cos(phi);
s1=sin(psi);
s2=sin(theta);
s3=sin(phi);
B = 1/cos(theta)*[0 s3 c3;
                  0 c2*c3 -c2*s3;
                  c2 s2*s3 s2*c3];

angles = B*omega;

dpsi = angles(1);
dtheta = angles(2);
dphi = angles(3);
end