function [T_stand,T_map_basic,T_Saas,T_Hop] = tropos(el,zd,t)
phi = 39.9929;
hd = 43000;
hw = 12000;
P0 = [1015.6 1016.6 1016.3 1018.6];
T0 = [283.7 289.8 298.15 294.5];
RH = [0.86 0.6 0.24 0.29];
e0 = 6.108.*RH*exp(1);
m = 1/(sqrt(1 - (cosd(el)/1.001)^2));
T_stand = m*zd;
[md,mw] = calcm(el);
T_map_basic = 2*md + 0.1*mw;
if t < 194400
    T_SaasD = 0.002277*(1+0.0023*cosd(2*phi)+0.00028*1.6485)*P0(1);
    T_SaasW = 0.002277*((1255/T0(1))+0.05)*e0(1);
    T_Saas = T_SaasD*md + T_SaasW*mw;
    T_HopD = 77.6e-6*(P0(1)/T0(1))*(hd/5);
    T_HopW = 0.373*(e0(1)/T0(1)^2)*(hw/5);
    T_Hop = T_HopD*md + T_HopW*mw;
elseif t < 216000
    T_SaasD = 0.002277*(1+0.0023*cosd(2*phi)+0.00028*1.6485)*P0(2);
    T_SaasW = 0.002277*((1255/T0(2))+0.05)*e0(2);
    T_Saas = T_SaasD*md + T_SaasW*mw;
    T_HopD = 77.6e-6*(P0(2)/T0(2))*(hd/5);
    T_HopW = 0.373*(e0(2)/T0(2)^2)*(hw/5);
    T_Hop = T_HopD*md + T_HopW*mw;
elseif t < 237600
    T_SaasD = 0.002277*(1+0.0023*cosd(2*phi)+0.00028*1.6485)*P0(3);
    T_SaasW = 0.002277*((1255/T0(3))+0.05)*e0(3);
    T_Saas = T_SaasD*md + T_SaasW*mw;
    T_HopD = 77.6e-6*(P0(3)/T0(3))*(hd/5);
    T_HopW = 0.373*(e0(3)/T0(3)^2)*(hw/5);
    T_Hop = T_HopD*md + T_HopW*mw;
else
    T_SaasD = 0.002277*(1+0.0023*cosd(2*phi)+0.00028*1.6485)*P0(4);
    T_SaasW = 0.002277*((1255/T0(4))+0.05)*e0(4);
    T_Saas = T_SaasD*md + T_SaasW*mw;
    T_HopD = 77.6e-6*(P0(4)/T0(3))*(hd/5);
    T_HopW = 0.373*(e0(4)/T0(4)^2)*(hw/5);
    T_Hop = T_HopD*md + T_HopW*mw;
end