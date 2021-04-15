function dxdt = EKF_ode(t,Z,n,w)
eta = 1000;
x = Z(1);
xdot = Z(2);
Phi_flat = Z(n+1:(n^2 + n));
Phi = reshape(Phi_flat,n,n);

Amat = [0 1;(-1 - 3*eta*x^2) 0];

Phi_dot = Amat*Phi;
Phi_dot_flat = reshape(Phi_dot,n^2,1);

dxdt = [xdot;(-1*x - eta*x^3 + w);Phi_dot_flat];
% dxdt = reshape(dxdt,2*N,1);
end