function yout = odefun(t,y)

h_0 = y(1);
v_0 = y(2);
m_0 = y(3);

if t < .26
    f_thrust = 60*t;
    dmdt = -0.01818;
elseif t < 1.65
    f_thrust = 15;
    dmdt = -.01818;
elseif t >= 1.65
    dmdt = 0;
    f_thrust = 0;
end

dvdt = (f_thrust / m_0) - 9.81 - (.000526 / m_0) * sign(v_0)*(v_0)^2;

dhdt = v_0;

yout = transpose([dhdt dvdt dmdt]);
        
    