function Q_DT = calc_DT_PN(dt,sigma)
sigmax = sigma(1);
sigmay = sigma(2);
sigmaz = sigma(3);

Q_DT(1,:) = [(dt^3*sigmax^2)/3, 0, 0, (dt^2*sigmax^2)/2, 0, 0, 0];
Q_DT(2,:) = [0, (dt^3*sigmay^2)/3, 0, 0, (dt^2*sigmay^2)/2, 0, 0];
Q_DT(3,:) = [0, 0, (dt^3*sigmaz^2)/3, 0, 0, (dt^2*sigmaz^2)/2, 0];
Q_DT(4,:) = [(dt^2*sigmax^2)/2, 0, 0, dt*sigmax^2, 0, 0, 0];
Q_DT(5,:) = [0, (dt^2*sigmay^2)/2, 0, 0, dt*sigmay^2, 0, 0];
Q_DT(6,:) = [0, 0, (dt^2*sigmaz^2)/2, 0, 0, dt*sigmaz^2, 0];
Q_DT(7,:) = [0, 0, 0, 0, 0, 0, 0];
end