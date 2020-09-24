function [out] = funcgrab2(j,t,u1,u2)
if j == 1
    out = (1/9)*u1 - (2/3)*u2 - (1/9)*t^2 + (2/3);
end
if j == 2
    out = u2 + 3*t - 4;
end
end