function dxdt = lindyn_ode(t,state,k,N,w)
state = reshape(state,N,2);
x = state(:,1);
xdot = state(:,2);

dxdt = [xdot;-k*x];
% dxdt = reshape(dxdt,2*N,1);
end