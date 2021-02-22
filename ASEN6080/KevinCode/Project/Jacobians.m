% ASEN 6080 - Project I
% Kevin Vose

%=================== symbolic differentiation =========================
syms x y z xdot ydot zdot xddot yddot zddot r mu Re J2 J3 omegaE
syms xs ys zs rho CD A m xs1 ys1 zs1 xs2 ys2 zs2 xs3 ys3 zs3
syms r0 H rho0 vxatm vyatm vzatm

% measurement partials
% setup s/c and station states
R = [x y z];
Rs = [xs ys zs];
V = [xdot ydot zdot];
dxs = -ys*omegaE;
dys = xs*omegaE;
dzs = 0;
Vs = [dxs dys dzs];
dxs1 = -ys1*omegaE;
dys1 = xs1*omegaE;
dzs1 = 0;
dxs2 = -ys2*omegaE;
dys2 = xs2*omegaE;
dzs2 = 0;
dxs3 = -ys3*omegaE;
dys3 = xs3*omegaE;
dzs3 = 0;

% measurements
range = norm(R - Rs);
dR = R - Rs;
dV = V - Vs;
drho = (dR(1)*dV(1)+dR(2)*dV(2)+dR(3)*dV(3))/range;

% partials w.r.t. x_sc (s/c state)
Y = [range; drho];

% partials w.r.t. x_obs (observer station state)
Robs = [xs ys zs];
H_obs_sc = jacobian(Y,Robs);

% ===============================================================
r = (x^2 + y^2 + z^2)^(1/2);   % radius of satellite orbit

Vatm = cross([0;0;omegaE],[x;y;z]);	% atmosphere velocity, ECI
vxatm = Vatm(1);
vyatm = Vatm(2);
vzatm = Vatm(3);

% velocity of s/c relative to atmosphere (in ECI)
Vrel = ((xdot - vxatm)^2 + (ydot - vyatm)^2 + (zdot - vzatm)^2)^(1/2);

% density of atmosphere
rho = rho0*exp(-(r-r0)/H);

% point mass accelerations
xddot_u = -mu*x/r^3;
yddot_u = -mu*y/r^3;
zddot_u = -mu*z/r^3;

% J2 accelerations
xddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(z/r)^2)*x/r);
yddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((1 - 5*(z/r)^2)*y/r);
zddot_J2 = -3/2*J2*(mu/r^2)*(Re/r)^2*((3 - 5*(z/r)^2)*z/r);

% drag accelerations
xddot_FD = -.5*rho*(CD*A/m)*Vrel*(xdot - vxatm);
yddot_FD = -.5*rho*(CD*A/m)*Vrel*(ydot - vyatm);
zddot_FD = -.5*rho*(CD*A/m)*Vrel*(zdot - vzatm);

% equations of motion
xddot = xddot_u + xddot_J2 + xddot_FD;
yddot = yddot_u + yddot_J2 + yddot_FD;
zddot = zddot_u + zddot_J2 + zddot_FD;

% state and differentiated state (18x1)
X = [x y z xdot ydot zdot mu J2 CD xs1 ys1 zs1 xs2 ys2 zs2 xs3 ys3 zs3];
%Xdot = [xdot ydot zdot xddot yddot zddot 0 0 0 0 0 0 0 0 0 0 0 0];
Xdot = [xdot ydot zdot xddot yddot zddot 0 0 0 dxs1 dys1 dzs1 dxs2 dys2 dzs2 dxs3 dys3 dzs3];

% A matrix jacobian
Amat = jacobian(Xdot,X);

% Htilde
Htilde_sc = jacobian(Y,X(1:9)); % compute this part of Htilde every time
                                % append zeros and H_obs_sc where appropriate (depending on station ID)


