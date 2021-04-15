function Q_DT = calc_DT_PN_1DOF(dt,sigma)
Q_DT(1,:) = [(dt^3*sigma^2)/3, (dt^2*sigma^2)/2];
Q_DT(2,:) = [(dt^2*sigma^2)/2,       dt*sigma^2];
end