function [value, isterminal, direction] = myEvent(time,outputs)

value = outputs - 0.002;
isterminal = 1;   %stops integration
direction = 0;
end