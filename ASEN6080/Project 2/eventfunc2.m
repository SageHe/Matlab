function [value, isterminal, direction] = eventfunc2(t,Z,BN)
r = Z(1:3);
Shat = BN(1,:);
value = dot(r,Shat);
isterminal = 1;
direction = 0;
end