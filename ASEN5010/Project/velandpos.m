function [r,rdot] = velandpos(radius,omega,i,theta)
global thetadotLMO thetadotGMO
if radius < 4000
    thetadot = thetadotLMO;
else
    thetadot = thetadotGMO;
end
q = [omega,i,theta];
omega = q(1);
i = q(2);
theta = q(3);
ir = [radius 0 0]';
rdot = [0 thetadot*radius 0]';
C = Euler3132C([omega,i,theta]');
r = C'*ir;
rdot = C'*rdot;
end