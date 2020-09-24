%HW4 problem 2
clear;close all;clc
%Define Inertia matrix
I = [125 0 0;0 100 0;0 0 75];
%Define starting position (using Euler angles)
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
N = 500;
t = linspace(0,500,N);
dt = t(2) - t(1);
%Calcualte initial angular momentum H
w = [1 0 0]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
%Case 2
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
w = [0 1 0]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
%Case 3
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
w = [0 0 1]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
%Case 4
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
w = [1 0.1 0]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
%Case 5
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
w = [0.1 1 0]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
%Case 6
Eulers = [10 10 10]';
Eulers = deg2rad(Eulers);
w = [0.1 0 1]'; %Initial spin rate
H1 = w(1)*I(1,1)/sin(Eulers(2));
H2 = -w(2)*I(2,2)/(sin(Eulers(3))*cos(Eulers(3)));
H3 = -w(3)*I(3,3)/(cos(Eulers(3))*cos(Eulers(3)));
for i = 1:N
    psidot(i) = -H1*((sin(Eulers(3,i))^2/I(2,2)) + (cos(Eulers(3,i))^2/I(3,3)));
    thetadot(i) = (H2/2)*((1/I(3,3)) - (1/I(2,2)))*sin(2*Eulers(3,i))*cos(Eulers(2,i));
    phidot(i) = H3*((1/I(1,1)) - ((sin(Eulers(3,i))^2)/I(2,2)) - ((cos(Eulers(3,i))^2)/I(3,3)))*sin(Eulers(2,i));
    
    Eulers(1,i+1) = Eulers(1,i) + psidot(i)*dt;
    Eulers(2,i+1) = Eulers(2,i) + thetadot(i)*dt;
    Eulers(3,i+1) = Eulers(3,i) + phidot(i)*dt;
end
figure
subplot(3,1,1)
plot(Eulers(1,:))
xlabel('Timesteps')
ylabel('\psi')
title('Angle vs Time')
subplot(3,1,2)
plot(Eulers(2,:))
xlabel('Timesteps')
ylabel('\theta')
subplot(3,1,3)
plot(Eulers(3,:))
xlabel('Timesteps')
ylabel('\phi')
figure
subplot(3,1,1)
plot(psidot)
xlabel('Timesteps')
ylabel('$\dot{\psi}$ (rad/s)', 'Interpreter','latex')
title('Angular Rate vs Time')
subplot(3,1,2)
plot(thetadot)
xlabel('Timesteps')
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
subplot(3,1,3)
plot(phidot)
xlabel('Timesteps')
ylabel('$\dot{\phi}$ (rad/s)', 'Interpreter','latex')
