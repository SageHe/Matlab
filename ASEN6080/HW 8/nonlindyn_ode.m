function dxdt = nonlindyn_ode(t,state,k,N,w)
eta = 1000;
state = reshape(state,N,2);
x = state(:,1);
xdot = state(:,2);

dxdt = [xdot;-k*x - eta*x.^3 + w];
% dxdt = reshape(dxdt,2*N,1);
end