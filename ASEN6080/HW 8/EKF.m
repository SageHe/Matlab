%% HW 8 -- Particle filter. Build linear KF and EKF (very basic) to compare to performance of the Particle filter
% Problem 1 -- linear dynamcis, linear measurement -- Implement and run EKF
% with SNC, generate requested plots
clear all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
rng('default'); %fix random number seed for repeatability
load('EKF_Noisy_Meas')
load('truth_state_nonlin')
truth_state = x(:,1:2);
clear x
t = linspace(0,2*pi,50);
k = 1;
x0 = [0;1];
x_minus(:,1) = x0;
x_plus(:,1) = x0;
dt = t(2) - t(1);
n = 2;
w = 0.1^2*randn(1,numel(t));

P0 = diag([0.2^2 0.2^2]);
Ri = 0.1^2;
sigma = 0.1;

Q_DT = calc_DT_PN_1DOF(dt,sigma);

%initialize filter
P_plus(:,:,1) = P0;
x_plus(:,1) = x0;
for i = 2:numel(t)
    %read observation
    yi = y(i);
    %propagate
    Phi = eye(2);
    Phi_flat = reshape(Phi,n^2,1);
    Z = [x_plus(:,i-1);Phi_flat];
    tspan = [0 dt];
    [~,x] = ode45(@(t,Z) EKF_ode(t,Z,n,w(i)),tspan,Z,opts);
    x_minus(:,i) = x(end,1:2)';
    Phi = reshape(x(end,(n+1:end)),n,n);
    %time update
    P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi';% + Q_DT;
    %process observations
    r(i) = y(i) - x_minus(1,i);
    Hi = [1 0];
    K(:,i) = P_minus(:,:,i)*Hi'*pinv(Hi*P_minus(:,:,i)*Hi' + Ri);
    %measurement update
    x_plus(:,i) = x_minus(:,i) + K(:,i)*r(i);
    P_plus(:,:,i) = (eye(2) - K(:,i)*Hi)*P_minus(:,:,i)*(eye(2) - K(:,i)*Hi)' + K(:,i)*Ri*K(:,i)';
    covbound(i,:) = [2*sqrt(P_plus(1,1,i)) 2*sqrt(P_plus(2,2,i))];
end

x_plus = x_plus';

figure
subplot(1,2,1)
hold on
plot(t,(x_plus(:,1) - truth_state(:,1)))
plot(t,covbound(:,1),'r')
plot(t,-covbound(:,1),'r')
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($x$)','Interpreter','latex')
legend('State Error','2\sigma Covariance')
subplot(1,2,2)
hold on
plot(t,(x_plus(:,2) - truth_state(:,2)))
plot(t,covbound(:,2),'r')
plot(t,-covbound(:,2),'r')
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($\dot{x})$','Interpreter','latex')
legend('State Error','2\sigma Covariance')
sgtitle('State Estimate Errors VS Time, EKF')


figure
hold on
plot(t,truth_state(:,1))
plot(t,y)
plot(t,x_plus(:,1))
grid on
grid minor
xlabel('Time')
ylabel('Measurement ($x$)','Interpreter','latex')
legend('True Measurement','Noisy Measurement','Predicted Measurement')
title('True VS Noisy VS Predicted Measurement Value, KF')

figure
hold on
plot(t,r)
plot(t,covbound(:,1),'r')
plot(t,-covbound(:,1),'r')
grid on
grid minor
xlabel('Time')
ylabel('Measurement Residual')
legend('Measurement Residual','2\sigma Covariance')
title('Measurement Residual VS Time, EKF')


x_plus = x_plus';
t = linspace(0,2*pi,200);
for i = 1:numel(t)
    [V,D] = eig(2*sqrtm(P_plus(:,:,end)));
    ell(:,i) = V*[sqrt(D(1,1))*cos(t(i));sqrt(D(2,2))*sin(t(i))];
end
% plot(ell(1,:)+x_plus(1,end),ell(2,:)+x_plus(2,end))
x_plus = x_plus';

figure
hold on
plot(truth_state(:,1),truth_state(:,2),'*')
plot(x_plus(:,1),x_plus(:,2))
plot(x_plus(1,1),x_plus(1,2),'og','MarkerSize',10)
plot(x_plus(end,1),x_plus(end,2),'^m','MarkerSize',10)
plot(ell(1,:)+x_plus(end,1),ell(2,:)+x_plus(end,2))
axis equal
grid on
grid minor
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
legend('True Traj.','Estimated Traj.','Initial Guess','Final Solution','2\sigma Covariance Ellipse')
title('EKF Phase Space Trajectory Plot')
