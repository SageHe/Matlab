function dydt = propmars(t,X_mars)
dydt = zeros(6,1);

mu = 1.32712428e11;

Rm = X_mars(1:3);
Vm = X_mars(4:6);

rd = Vm;
rdd = (-mu.*Rm)./(norm(Rm))^3;

dydt(1:3) = rd;
dydt(4:6) = rdd;
end