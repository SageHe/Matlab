syms dt sigma w1 w2 w3

A = zeros(9,9);
A(1:6,4:end) = eye(6);
stm = eye(9) + A*dt + .5*A*A*dt^2;
B = [zeros(6,3);eye(3)];

assumeAlso(dt,'real')

Q = diag([w1 w2 w3]);

Q_DT = stm*B*Q*B'*stm';

Q_DT = int(Q_DT,dt);
