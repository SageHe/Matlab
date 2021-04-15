clear all;close all;clc
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
% HW 8 -- PF with importance sampling using initial guess for x and xdot
% from U(-2,2) and N = 5000
rng('default')
% load('Noisy_Meas')
load('truth_state_nonlin')
truth_state = x(:,1:2);
clear x;
%given info
a = -2;
b = 2;
kdyn = 1;
N = 5000;
t = linspace(0,2*pi,50);
dt = t(2) - t(1);
eta = 0.1^2*randn(1,numel(t));
y = truth_state(:,1) + eta';

% A = [0 1;-1 0];

%initialize filter
x0 = a + (b-a).*rand(N,2);
x = a + (b-a).*rand(N,2);
w(:,1) = 1/N*ones(N,1);
Ri = 0.1^2;
w_pn = normrnd(0,0.1,N,1);
sigma = 0.1^2;
Q_DT = calc_DT_PN_1DOF(dt,sigma);
% figure
% hold on
% plot(truth_state(:,1),truth_state(:,2))
for i = 2:numel(t)
   %obtain measurement 
   ybar_i = y(i);
   %propagate
%    for k = 1:N
%        A = [0 1;-1 0];
%        Phi = expm(A*dt);
%        x(k,:,i) = Phi*x(k,:,i-1)'; %+ sqrtm(Q_DT)*rand(2,1);
%        %measurement update
%        yhat(k,i) = x(k,1,i);
%        r(k,i) = ybar_i - yhat(k,i);
%        w(k,i) = 1/(((2*pi)^(1/2))*Ri^(1/2))*exp(-.5*r(k,i)*inv(Ri)*r(k,i));
%    end
%    x(:,:,i) = x(:,:,i) + mvnrnd([0 0], 0.01*[(1/3)*dt^3 (1/2)*dt^2; (1/2)*dt^2 dt]);
    %propagate
    tspan = [0 dt];
    x = reshape(x,2*N,1);
    [~,x] = ode45(@(t,x) nonlindyn_ode(t,x,kdyn,N,w_pn),tspan,x,opts);
    x = reshape(x(end,:),N,2);
%     x = x + mvnrnd([0 0], 0.01*[(1/3)*dt^3 (1/2)*dt^2; (1/2)*dt^2 dt],N);
    %measurement update
    for k = 1:N
        yhat(k,i) = x(k,1);
        r(k,i) = ybar_i - yhat(k,i);
        w(k,i) = 1/(((2*pi)^(1/2))*Ri^(1/2))*exp(-.5*r(k,i)*inv(Ri)*r(k,i));
    end
   %resample
   wsum = w(1,i);
   for k = 2:N
       wsum(k) = wsum(k-1) + w(k,i);
       wcdf = wsum./wsum(end);
   end
   for k = 1:N
       clear ind
       u(k) = rand(1);
       ind = find(wcdf < u(k));
       if u(k) < min(wcdf)
           j = 1;
       else
           j = ind(end) + 1;
       end
       xstar_k(k,:,i) = x(j,:);
       w(k,i) = 1/N;
   end
   x = xstar_k(:,:,i);
   xstar(i,:) = mean(xstar_k(:,:,i));
   P(:,:,i) = cov(xstar_k(:,:,i));
   covbound(i,:) = [2*sqrt(P(1,1,i)) 2*sqrt(P(2,2,i))];
   rplot(i) = ybar_i - xstar(i,1);
%    plot(xstar(:,1),xstar(:,2))
end
%% 
figure
subplot(1,2,1)
hold on
plot(t,(xstar(:,1) - truth_state(:,1)))
plot(t,covbound(:,1),'r')
plot(t,-covbound(:,1),'r')
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($x$)','Interpreter','latex')
legend('State Error','2\sigma Covariance')
subplot(1,2,2)
hold on
plot(t,(xstar(:,2) - truth_state(:,2)))
plot(t,covbound(:,2),'r')
plot(t,-covbound(:,2),'r')
ylim([-2 2])
grid on
grid minor
xlabel('Time')
ylabel('State Estimate Error ($\dot{x})$','Interpreter','latex')
legend('State Error','2\sigma Covariance')
sgtitle('State Estimate Errors VS Time, PF')

figure
hold on
plot(t,rplot)
plot(t,covbound(:,1),'r')
plot(t,-covbound(:,1),'r')
grid on
grid minor
xlabel('Time')
ylabel('Measurement Residual')
legend('Measurement Residual','2\sigma Covariance')
title('Measurement Residual VS Time, PF')

t = linspace(0,2*pi,200);
for i = 1:numel(t)
    [V,D] = eig(2*sqrtm(P(:,:,end)));
    ell(:,i) = V*[sqrt(D(1,1))*cos(t(i));sqrt(D(2,2))*sin(t(i))];
end
% plot(ell(1,:)+x_plus(1,end),ell(2,:)+x_plus(2,end))

figure
hold on
plot(truth_state(:,1),truth_state(:,2),'*')
plot(xstar(:,1),xstar(:,2))
plot(xstar(1,1),xstar(1,2),'og','MarkerSize',10)
plot(xstar(end,1),xstar(end,2),'^m','MarkerSize',10)
plot(ell(1,:)+xstar(end,1),ell(2,:)+xstar(end,2))
plot(xstar_k(:,1,end),xstar_k(:,2,end),'.b')
legend('True Traj.','Estimated Traj.','Initial Guess','Final Solution','2\sigma Covariance Ellipse','Particles @ t_f')
h = get(gca,'Children');
set(gca,'Children',[h(6) h(5) h(4) h(3) h(2) h(1)]);
grid on
grid minor
axis equal
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
% legend('True Traj.','Estimated Traj.','Initial Guess','Final Solution','2\sigma Covariance Ellipse','Particles @ t_f')
title('PF Phase Space Trajectory Plot')

% plot(xstar_k(:,1,i),xstar_k(:,2,i),'k.')
% t = linspace(0,2*pi,200);
% for i = 1:numel(t)
%     [V,D] = eig(sqrtm(P(:,:,end)));
%     ell(:,i) = V*[sqrt(D(1,1))*cos(t(i));sqrt(D(2,2))*sin(t(i))];
% end
% plot(ell(1,:),ell(2,:))
