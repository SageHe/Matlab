function [child] = imagebreed(par1,par2)
[a b c] = size(par1);
vec = [];
for i = 1:numel(par1)
    d = randi([0,1]);
    vec = [vec d];
    %vec(i) = a;
end
vec = reshape(vec,a,b,c);
child = par1 .* vec;
x = find(~child);
for j = 1:length(x)
    child(x(j)) = par2(x(j));
end
end
    
    