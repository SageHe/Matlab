function [out] = funcgrab4(j,t,u1,u2,u3)
if j == 1
    out = u2;
end
if j == 2
    out = u3;
end
if j == 3
    out = (8*t^3 - 2 - 2*u1 + 2*t*u2 - t^2*u3)/t^3;
end
end