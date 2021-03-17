
syms dt 

A = zeros(6,6);
A(1:3,4:6) = eye(3);

stm = eye(6) + A*dt;
B = [zeros(3);eye(3)];

assumeAlso(dt,'real')

intvar = stm*B; 

Gamma = int(intvar,dt);
