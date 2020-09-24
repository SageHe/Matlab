function [out] = prodofb(a,b)
if 2 * b - 1 > a
    out = 0;
else
    out = 1;
end
for i = 1:2:(2 * b - 1)
    out = out * i;
end
    