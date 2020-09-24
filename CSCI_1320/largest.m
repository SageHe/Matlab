function out = largest(vec)
if max(vec) <= 0
    out = [];
    return;
else if isempty(vec) == 1
    out = [];
    return;
else if isequal(find(vec == 0), [])
    out = [];
    return;
else if sum(vec) == 0
    out = [];
    return;
    end
    end
    end
end
zero = find(vec == 0);
right_shift = zero + 1;
left_shift = zero - 1;
if right_shift(end) > length(vec)
    right_shift(end) = [];
end
if left_shift(1) < 1
    left_shift(1) = [];
end
tot = [left_shift right_shift];
num = vec(tot);
out = max(num);
end

