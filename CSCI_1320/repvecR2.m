function out = repvecR2(vec,x)
a = 1;
if x > 1
mat = [vec;vec];
end 
if x == 1
    mat = vec;
end
if x == 0
    mat = [];
end
while numel(mat) < (numel(vec) * x)
    mat = [mat;vec];
end
out = reshape(mat,1,numel(mat));
end