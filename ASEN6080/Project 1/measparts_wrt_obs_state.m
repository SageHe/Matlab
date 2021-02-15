syms x y z xdot ydot zdot xs ys zs xsdot ysdot zsdot

R = [x y z];
Rs = [xs ys zs];
V = [xdot ydot zdot];
Vs = [xsdot ysdot zsdot];

rho = norm(R - Rs);
rdiff = R - Rs;
vdiff = V - Vs;
num = (rdiff(1)*vdiff(1)) + (rdiff(2)*vdiff(2)) + (rdiff(3)*vdiff(3));

rhodot = num/rho;

obs_state = [Rs];

f = [rho rhodot];

rdd = simplify(jacobian(f,obs_state));

