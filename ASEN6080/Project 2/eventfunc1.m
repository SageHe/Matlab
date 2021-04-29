function [value, isterminal, direction] = eventfunc1(t,Z,opts)
RSOI = 925000;
value = norm(Z(1:3)) - 5*RSOI;
isterminal = 1;
direction = 0;
end