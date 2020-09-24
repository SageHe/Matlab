function out = causeMutation(child)
%tic
global mutation_rate
ascii = [32 65:90 97:122];
%mutation_rate = input('What is the mutation rate(percentage)?');
vec = 1:(round(100/mutation_rate));%If mutation rate is input as 20 percent, vec becomes a 5 element vector
for i = 1:length(child)%Counter from 1 to length of child
    if vec(randi(length(vec))) == 1%If the randomly chosen element from vec is a one, mutate the child
        child(i) = ascii(randi(length(ascii)));
    end
end
out = child;
%toc
end
    