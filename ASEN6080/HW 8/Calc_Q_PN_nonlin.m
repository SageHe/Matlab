
syms dt sigma

A = [0 1;-1 0];

stm = eye(2) + A*dt;
B = [0;1];

assumeAlso(dt,'real')

Q = sigma^2;

Q_DT = stm*B*Q*B'*stm';

Q_DT = int(Q_DT,dt);
