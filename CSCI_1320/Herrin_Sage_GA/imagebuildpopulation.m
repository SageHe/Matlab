function [population] = imagebuildpopulation(target,population_size)
    [a b c] = size(target);
    population = randi(256,a,b,c,population_size); %create a 4D matrix, with n 'layers' of test images all the same sizes as the target
end