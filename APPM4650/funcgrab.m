function [out] = funcgrab(j,t,u1,u2)
if j == 1
    out = u1 - u2 + 2;
end
if j == 2
    out = -u1 + u2 + 4*t;
end
end