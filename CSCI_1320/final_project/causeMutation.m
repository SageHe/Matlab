function out = causeMutation(child)
%tic
global mutation_rate
ascii = [32 65:90 97:122];
%mutation_rate = input('What is the mutation rate(percentage)?');
vec = 1:(round(100/mutation_rate));
for i = 1:length(child)
    if vec(randi(length(vec))) == 1
        child(i) = ascii(randi(length(ascii)));
    end
end
out = child;
%toc
end
    