function Q_DT = calc_DT_PN_DMC(dt,wvec)
w1 = wvec(1);
w2 = wvec(2);
w3 = wvec(3);

Q_DT(1,:) = [(dt^5*w1)/20, 0, 0, (dt^4*w1)/8, 0, 0, (dt^3*w1)/6, 0, 0];
Q_DT(2,:) = [0, (dt^5*w2)/20, 0, 0, (dt^4*w2)/8, 0, 0, (dt^3*w2)/6, 0];
Q_DT(3,:) = [0, 0, (dt^5*w3)/20, 0, 0, (dt^4*w3)/8, 0, 0, (dt^3*w3)/6];
Q_DT(4,:) = [(dt^4*w1)/8, 0, 0, (dt^3*w1)/3, 0, 0, (dt^2*w1)/2, 0, 0];
Q_DT(5,:) = [0, (dt^4*w2)/8, 0, 0, (dt^3*w2)/3, 0, 0, (dt^2*w2)/2, 0];
Q_DT(6,:) = [0, 0, (dt^4*w3)/8, 0, 0, (dt^3*w3)/3, 0, 0, (dt^2*w3)/2];
Q_DT(7,:) = [(dt^3*w1)/6, 0, 0, (dt^2*w1)/2, 0, 0, dt*w1, 0, 0];
Q_DT(8,:) = [0, (dt^3*w2)/6, 0, 0, (dt^2*w2)/2, 0, 0, dt*w2, 0];
Q_DT(9,:) = [0, 0, (dt^3*w3)/6, 0, 0, (dt^2*w3)/2, 0, 0, dt*w3];
end