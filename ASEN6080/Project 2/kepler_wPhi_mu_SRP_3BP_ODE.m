function zd = kepler_wPhi_mu_SRP_3BP_ODE(t,Z,const)

%Setting up function params
n = const.n;
mu = const.mu;
mu_s = const.mu_s;
Am = const.Am;
P_Phi = const.P_Phi;
% Phi_SRP = const.Phi_SRP;
JD_o = const.JD_o;
AU = const.AU;
% c = const.c;

% Getting Earth Orbit info
JD = JD_o + t/86400;
% RsE = planetEphemeris(JD,'earth','sun')';
T = (JD - 2451545.0)/36525;
deg2rad = pi/180;
% Earth ephemeris coefficients
coeff.L = [100.466449     35999.3728519  -0.00000568      0.0];
coeff.a = [1.000001018    0.0             0.0             0.0];
coeff.e = [0.01670862    -0.000042037    -0.0000001236	0.00000000004];
coeff.i = [0.0            0.0130546      -0.00000931     -0.000000034];
coeff.W = [174.873174    -0.2410908       0.00004067     -0.000001327];
coeff.P = [102.937348     0.3225557       0.00015026      0.000000478];
coeff.mu_p = 3.98600432896939e5;
Tvec = [1; T; T*T; T*T*T];
% Mean longitude of Planet
L = coeff.L*Tvec*deg2rad;
% Semimajor axss of the orbit
a = coeff.a*Tvec*AU;
% Eccentricity of the orbit
e = coeff.e*Tvec;
% Inclination of the orbit
inc = coeff.i*Tvec*deg2rad;
% Longitude of the Ascending Node
W = coeff.W*Tvec*deg2rad;
% Longitude of the Perihelion
P = coeff.P*Tvec*deg2rad;
% Argument of perihelion
w = P - W;
% Mean anomaly of orbit
M = L - P;
% True anomaly of orbit
Ccen = (2*e - e^3/4 + 5/96*e^5)*sin(M) + (5/4*e^2 - 11/24*e^4)*sin(2*M) + ...
    (13/12*e^3 - 43/64*e^5)*sin(3*M) + 103/96*e^4*sin(4*M) + ...
    1097/960*e^5*sin(5*M);
nu = M + Ccen;
% Trajectory equation
r = a*(1-e^2)/(1+e*cos(nu));
% Obtain position in sun orbital frame 
o = r*[cos(nu) sin(nu) 0];
% Transform Earth rel. to Sun position to inertial frame
xs = -(o(1)*(cos(w)*cos(W)-sin(w)*cos(inc)*sin(W))-o(2)*(sin(w)*cos(W)+cos(w)*cos(inc)*sin(W)));
ys = -(o(1)*(cos(w)*sin(W)+sin(w)*cos(inc)*cos(W))+o(2)*(cos(w)*cos(inc)*cos(W)-sin(w)*sin(W)));
zs = -(o(1)*sin(w)*sin(inc)+o(2)*cos(w)*sin(inc));
Rs = [xs;ys;zs];
% Convert to EME2000
theta = 23.4393*deg2rad;
    C = [1 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)];
    Rs = C*Rs;
xs = Rs(1); ys = Rs(2); zs = Rs(3);

% Unpack State Vector
x = Z(1); y = Z(2); z = Z(3);
v_vec = Z(4:6); 
Cr = Z(7);
Phi_flat = Z(n+1:end);
Phi = reshape(Phi_flat,n,n);

% Compute Derivatives
a_vec = [ 
(AU^2*Am*Cr*P_Phi*(x - xs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (mu*x)/(x^2 + y^2 + z^2)^(3/2) - mu_s*(xs/(xs^2 + ys^2 + zs^2)^(3/2) + (x - xs)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2));
(AU^2*Am*Cr*P_Phi*(y - ys))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (mu*y)/(x^2 + y^2 + z^2)^(3/2) - mu_s*((y - ys)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) + ys/(xs^2 + ys^2 + zs^2)^(3/2));
(AU^2*Am*Cr*P_Phi*(z - zs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (mu*z)/(x^2 + y^2 + z^2)^(3/2) - mu_s*((z - zs)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) + zs/(xs^2 + ys^2 + zs^2)^(3/2))];

Xdot_vec = [v_vec; a_vec; 0];
A_t = [
                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                    0,                                                                                                                                                                                                                                                                                                                                                                    0, 1, 0, 0,                                                                     0;
                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                    0,                                                                                                                                                                                                                                                                                                                                                                    0, 0, 1, 0,                                                                     0;
                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                    0,                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 1,                                                                     0;
(3*mu*x^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*x - 2*xs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am*Cr*P_Phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am*Cr*P_Phi*(2*x - 2*xs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)),                                                                                                                                                   (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*y - 2*ys)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*y - 2*ys)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)),                                                                                                                                                   (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*z - 2*zs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*z - 2*zs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am*P_Phi*(x - xs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2);
                                                                                                                                                  (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*x - 2*xs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*x - 2*xs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*y^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*y - 2*ys)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am*Cr*P_Phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am*Cr*P_Phi*(2*y - 2*ys)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)),                                                                                                                                                   (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*z - 2*zs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*z - 2*zs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am*P_Phi*(y - ys))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2);
                                                                                                                                                  (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*x - 2*xs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*x - 2*xs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)),                                                                                                                                                   (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*y - 2*ys)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am*Cr*P_Phi*(2*y - 2*ys)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*z^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*z - 2*zs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am*Cr*P_Phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am*Cr*P_Phi*(2*z - 2*zs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am*P_Phi*(z - zs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2);
                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                    0,                                                                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0,                                                                     0];
Phi_dot = A_t*Phi;

%Reformat for Output
Phi_dot_flat = reshape(Phi_dot,n^2,1);
zd = [Xdot_vec; Phi_dot_flat];
end