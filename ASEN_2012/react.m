function dydt = react(t,y)

k1 = 5;
k2 = 2;
k3 = 1;

A = y(1);
B = y(2);
C = y(3);

dydt(1)=-k1*A+k2*B;
dydt(2)=k1*A-(k2+k3)*B;
dydt(3)=k3*B;

dydt = dydt';

end