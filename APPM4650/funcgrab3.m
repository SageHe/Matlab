function [out] = funcgrab3(j,t,u1,u2,u3)
if j == 1
    out = u2;
end
if j == 2
    out = u3;
end
if j == 3
    out = -u3 +4*u2 + 4*u1;
end
end