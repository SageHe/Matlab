function [newchild] = imagemutation(child, mutationrate, range, randomrate)
   [a b c] = size(child);%Creates variables a, b, and c according to size of child image
    mutation = randi(round((100/mutationrate)),a,b,c);%Depending on the input mutation rate(percentage), creates an array of values proportional size depending on mutation rate, randomly picks value from array, and if random value is one, mutates child
    x = find(mutation == 1);
   for i = 1:length(x)
       randommutation = randi(round(100/randomrate));%Creates another chance that if the child is radnomly mutated that that random mutation is completely random, not inside the mutation range
       if randommutation == 1
           child(x(i)) = randi(256);
       else
           child(x(i)) = child(x(i)) + randi([-range, range]);
       end
   end
   newchild = child;
end
    
       