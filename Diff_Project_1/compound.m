function [out] = compound(n,t,Y_0)
    y = [];
    for i = 1:length(t)
    y = [y (1 + (.03/n))^(n*t(i))*Y_0];
    end
    out = y;
end