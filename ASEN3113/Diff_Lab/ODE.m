tspan = [0 500];
y0 = [15 1 0.1];

%Call ode45 and get solution
[t,y] = ode45(@Rossler,tspan,y0);

plot3(y(:,1),y(:,2),y(:,3))
title('Rossler Attractor')
