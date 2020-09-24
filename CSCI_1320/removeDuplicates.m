function out = removeDuplicates(vec)
long = [];
for i = 1:max(vec)
    A = (find(vec == i));
    B = length(A);
    if B >= 1
        long = [long i];
        vec(A) = [];
    end
end
out = long;
out = length(long);
end
        
        
        

