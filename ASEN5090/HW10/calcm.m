function [md,mw] = calcm(el)
md = 1/(sind(el) + (0.00143/(tand(el) + 0.0445)));
mw = 1/(sind(el) + (0.00035/(tand(el) + 0.017)));
end