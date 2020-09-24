function [newchild] = imagemutation(child, mutationrate, range, randomrate)
   [a b c] = size(child);
    mutation = randi(round((100/mutationrate)),a,b,c);
    x = find(mutation == 1);
   for i = 1:length(x)
       randommutation = randi(round(100/randomrate));
       if randommutation == 1
           child(x(i)) = randi(256);
       else
           child(x(i)) = child(x(i)) + randi([-range, range]);
       end
   end
   newchild = child;
end
    
       