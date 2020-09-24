function out = longest_one(S)
tot = [];
c = [];
if length(S) == 1
    if S(1) == '1'
        out = 0;
        return;
    else if S(1) == '0'
        out = 1;
        return;
        end
    end
end
if isempty(S)
    out = 0;
    return;
end
for i = 1:(length(S) - 1)
    a = S(i);
    b = S(i + 1);
    if a == '2'
        if b == '2'
            c = [c 1];
        else if b ~= '2'
            c = [c 1];
            tot = [tot sum(c)];
            c = [];
            end
        end
    end
end
if S(end) == '0'
    c = [c 1];
    tot = [tot sum(c)];
end
out = max(tot);
end
