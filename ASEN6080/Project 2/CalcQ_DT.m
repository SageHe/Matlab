
syms dt sigma

A = zeros(7,7);
A(1:3,4:6) = eye(3);

stm = eye(7) + A*dt;
B = [zeros(3);eye(3);zeros(1,3)];

assumeAlso(dt,'real')

Q = sigma^2*eye(3);

Q_DT = stm*B*Q*B'*stm';

Q_DT = int(Q_DT,dt);
