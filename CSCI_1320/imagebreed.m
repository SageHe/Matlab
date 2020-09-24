function [child] = imagebreed(par1,par2)
[a b c] = size(par1);%Creates variable a, b, and c of corresponding size of the parent image
vec = [];
for i = 1:numel(par1)%Counter from one to the number of elements in the parent vector
    d = randi([0,1]);%Generates a random zero or a one
    vec = [vec d];%Concatenates vec and the randomly generated value of d
    %vec(i) = a;
end
vec = reshape(vec,a,b,c);%reshpaces vec into an array of the same dimensions as the parent
child = par1 .* vec;%Assigns random parts of the child the same elements as parent 1
x = find(~child);
for j = 1:length(x)%Assigns all other parts of the child to be the same as parent 2
    child(x(j)) = par2(x(j));
end
end
    
    