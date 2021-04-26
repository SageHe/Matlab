function dxdt = TBG_SRP_ode(t,state,n)
JD0 = 2456296.25;
c = 299792.458; %speed of light in km/s
mu_s = 132712440017.987; %grav. param of the sun in km^3/s^2
Am_rat = 0.01/(1e6); %area to mass ratio in km^2/kg
Cr = state(7); %Cannonball SRP coefficient
AU = 149597870.700; %AU in km

JD = JD0 + t/86400;

x = state(1); %position and velocity of s/c wrt to Earth
y = state(2);
z = state(3);
R_sc = [x y z]';
R_sc_norm = norm(R_sc);
xdot = state(4);
ydot = state(5);
zdot = state(6);
V_sc = [xdot ydot zdot]';

STM_flat = state(n+1:end);
STM = reshape(STM_flat,n,n);

[Re_s,Ve_s,mu_p] = Ephem(JD,3,'EME2000'); %gives position of Earth wrt to the Sun
mu = mu_p;
Rs_e = -Re_s; %negate vector in order to get position of sun wrt to Earth
xs = Rs_e(1); ys = Rs_e(2);zs = Rs_e(3);
Rs_e = [xs ys zs]';
Rs_e_norm = norm(Rs_e);
Vs_e = -Ve_s;
ri_sc = Rs_e - R_sc;
ri_sc_norm = norm(ri_sc);
uhat = ri_sc/ri_sc_norm;

Phi = 1357;
P_phi = Phi/c;
%calculate acceleration due to SRP
a_srp = -Cr*(P_phi/ri_sc_norm^2)*Am_rat*AU^2*uhat;
%calculate acceleration due to 3BG
a_3BG = mu_s*((ri_sc/ri_sc_norm^3) - (Rs_e/Rs_e_norm^3));
%calculate normal point mass acceleration due to Earth
a_kep(1,1) = -(mu*x)/(x^2 + y^2 + z^2)^(3/2);
a_kep(2,1) = -(mu*y)/(x^2 + y^2 + z^2)^(3/2);
a_kep(3,1) = -(mu*z)/(x^2 + y^2 + z^2)^(3/2);
%add accelerations together for total acceleration
a = a_kep + a_srp + a_3BG;

Amat(1,:) = [0, 0, 0, 1, 0, 0, 0];
Amat(2,:) = [0, 0, 0, 0, 1, 0, 0];
Amat(3,:) = [0, 0, 0, 0, 0, 1, 0];
Amat(4,:) = [(3*mu*x^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*x - 2*xs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am_rat*Cr*P_phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am_rat*Cr*P_phi*(2*x - 2*xs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*y - 2*ys)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*y - 2*ys)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*z - 2*zs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*z - 2*zs)*(x - xs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am_rat*P_phi*(x - xs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2)];
Amat(5,:) = [(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*x - 2*xs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*x - 2*xs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*y^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*y - 2*ys)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am_rat*Cr*P_phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am_rat*Cr*P_phi*(2*y - 2*ys)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*z - 2*zs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*z - 2*zs)*(y - ys))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am_rat*P_phi*(y - ys))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2)];
Amat(6,:) = [(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*x - 2*xs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*x - 2*xs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) + (3*mu_s*(2*y - 2*ys)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)) - (3*AU^2*Am_rat*Cr*P_phi*(2*y - 2*ys)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), (3*mu*z^2)/(x^2 + y^2 + z^2)^(5/2) - mu_s*(1/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*(2*z - 2*zs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2))) - mu/(x^2 + y^2 + z^2)^(3/2) + (AU^2*Am_rat*Cr*P_phi)/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2) - (3*AU^2*Am_rat*Cr*P_phi*(2*z - 2*zs)*(z - zs))/(2*((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(5/2)), 0, 0, 0, (AU^2*Am_rat*P_phi*(z - zs))/((x - xs)^2 + (y - ys)^2 + (z - zs)^2)^(3/2)];
Amat(7,:) = [0, 0, 0, 0, 0, 0, 0];

STM_dot = Amat*STM;
STM_dot_flat = reshape(STM_dot,n^2,1);

dxdt = [xdot;ydot;zdot;a;0;STM_dot_flat];
end
