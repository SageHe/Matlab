function [NwRcN] = calcNwRcN(t)
dt = 1e-3;
RcN331 = calcRcN(t+dt);
RcN330 = calcRcN(t);

RcNdot = (RcN331-RcN330)/dt;
omegatilde = -RcNdot*RcN330';
NwRcN = RcN330'*[mean([omegatilde(3,2),-omegatilde(2,3)]);
    mean([-omegatilde(3,1),omegatilde(1,3)]);
    mean([omegatilde(2,1),-omegatilde(1,2)])];
end