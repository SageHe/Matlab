%% HW 8 -- Particle filter. Build linear KF and EKF (very basic) to compare to performance of the Particle filter
% Problem 1 -- linear dynamcis, linear measurement -- Implement and run KF
% with SNC, generate requested plots
clear all;close all;clc
rng('default'); %fix random number seed for repeatability
t = linspace(0,2*pi,50);
k = 1;
A = [0 1;-1 0];
x0 = [0;1];
x(:,1) = x0;
dt = t(2) - t(1);

P0 = diag([0.2^2 0.2^2]);
Ri = 0.1^2;
sigma = 0.1^2;

P_plus(:,:,1) = P0;
eta = 0.1^2*randn(1,numel(t));
x_minus(:,1) = x0;
x_plus(:,1) = x0;
for i = 2:numel(t)
    Phi = expm(A*dt);
    x(:,i) = Phi*x(:,i-1);
    
    y(i) = x(1,i) + eta(i);
    %time update
    x_minus(:,i) = Phi*x_plus(:,i-1);
    
    Q_DT = calc_DT_PN_1DOF(dt,sigma);
    
    P_minus(:,:,i) = Phi*P_plus(:,:,i-1)*Phi' + Q_DT;
    %process observations
    r(i) = y(i) - x_minus(1,i);
    Hi = [1 0];
    K(:,i) = P_minus(:,:,i)*Hi'*pinv(Hi*P_minus(:,:,i)*Hi' + Ri);
    %measurement update
    x_plus(:,i) = x_minus(:,i) + K(:,i)*(y(i) - Hi*x_minus(:,i));
    P_plus(:,:,i) = (eye(2) - K(:,i)*Hi)*P_minus(:,:,i)*(eye(2) - K(:,i)*Hi)' + K(:,i)*Ri*K(:,i)';
    covbound(i,:) = [2*sqrt(P_plus(1,1,i)) 2*sqrt(P_plus(2,2,i))];
%     r(i) = y(i) - x_plus(1,i);
end

figure
subplot(1,2,1)
hold on
plot(t,(x_plus(1,:) - x(1,:)))
plot(t,covbound(:,1),'r')
plot(t,-covbound(:,1),'r')
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($x$)','Interpreter','latex')
legend('State Error','2\sigma Covariance')
subplot(1,2,2)
hold on
plot(t,(x_plus(2,:) - x(2,:)))
plot(t,covbound(:,2),'r')
plot(t,-covbound(:,2),'r')
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($\dot{x})$','Interpreter','latex')
legend('State Error','2\sigma Covariance')
sgtitle('State Estimate Errors VS Time, KF')

figure
hold on
plot(t,x(1,:))
plot(t,y)
plot(t,x_plus(1,:))
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
title('Measurement Residual VS Time, KF')

t = linspace(0,2*pi,200);
for i = 1:numel(t)
    [V,D] = eig(2*sqrtm(P_plus(:,:,end)));
    ell(:,i) = V*[sqrt(D(1,1))*cos(t(i));sqrt(D(2,2))*sin(t(i))];
end
% plot(ell(1,:)+x_plus(1,end),ell(2,:)+x_plus(2,end))

figure
hold on
plot(x(1,:),x(2,:),'*')
plot(x_plus(1,:),x_plus(2,:))
plot(x_plus(1,1),x_plus(2,1),'og','MarkerSize',10)
plot(x_plus(1,end),x_plus(2,end),'^m','MarkerSize',10)
plot(ell(1,:)+x_plus(1,end),ell(2,:)+x_plus(2,end))
grid on
grid minor
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
legend('True Traj.','Estimated Traj.','Initial Guess','Final Solution','2\sigma Covariance Ellipse')
title('KF Phase Space Trajectory Plot')