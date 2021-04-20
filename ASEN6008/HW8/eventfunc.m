function [value, isterminal, direction] = eventfunc(t,x,opts)
value = x(2);
isterminal = 1;
direction = 1;
end