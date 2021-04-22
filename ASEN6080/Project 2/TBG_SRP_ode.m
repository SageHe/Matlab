function dxdt = TBG_SRP_ode(t,state)
JD0 = 2456296.25;
c = 299792.458;
mu_S = 32712440017.987;

JD = JD0 + t/86400;

x = state(1);
y = state(2);
z = state(3);
xdot = state(4);
ydot = state(5);
zdot = state(6);



dxdt = [x*t;y;z;xdot*t;ydot;zdot;0];
end
