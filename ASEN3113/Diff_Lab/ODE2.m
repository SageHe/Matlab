tspan = [0,100];

y0 =[1;2;15];

[t,y] = ode45(@lorentz,tspan,y0);

hold on
plot(t,y(:,1))
plot(t,y(:,2));
hold off

figure
plot(y(:,1),y(:,2));