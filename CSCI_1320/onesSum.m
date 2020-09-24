function [out] = onesSum(a)
b = 0;
vec = [];
for i = 1:a
    b = b + (1 * 10^(i - 1));
    vec = [vec b];
end
out = sum(vec);
end
