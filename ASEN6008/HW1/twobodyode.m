function dydt = twobodyode(t,X)
dydt = zeros(36,1);

mu = 1.32712428e11;
mu_e = 3.986004415e5;
mu_m = 4.305e4;

R_sc = X(1:3);
V_sc = X(4:6);
Re = X(7:9);
Ve = X(10:12);
Rm = X(13:15);
Vm = X(16:18);

R_sc_e = R_sc - Re;
R_sc_m = R_sc - Rm;

rd_sc = V_sc;
rdd_sc = (-mu.*R_sc)./(norm(R_sc))^3;
rd_e = Ve;
rdd_e = (-mu.*Re)./(norm(Re))^3;
rd_m = Vm;
rdd_m = (-mu.*Rm)./(norm(Rm))^3;

dydt(1:3) = rd_sc;
dydt(4:6) = rdd_sc;
dydt(7:9) = rd_e;
dydt(10:12) = rdd_e;
dydt(13:15) = rd_m;
dydt(16:18) = rdd_m;

rd_sc_pert = V_sc;
rdd_sc_pert = (-mu.*R_sc)./(norm(R_sc))^3 + 
end