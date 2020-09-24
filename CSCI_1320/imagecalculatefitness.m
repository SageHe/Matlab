function [fitness_scores] = imagecalculatefitness(target,population,tolerance,exponential,population_size)
    [a b c] = size(target);%Creates variables a, b, and c of corresponding size to the size of the image
    pop_fitness = zeros(population_size,1);%Makes pop_fitness an array of all zeros
    for i = 1:population_size%Counter from one to the size of the population
        pop_mem = population(:,:,:,i);%Finds where the current population member is equal to the target image to within some tolerance(difference) and creates a logical vetor of where the population member is and is not fit
        pop_mem_fitness = sum(sum(sum((abs(target - pop_mem))<= tolerance)));
        pop_mem_fitness = pop_mem_fitness / numel(pop_mem);
        pop_fitness(i)= pop_mem_fitness;
        pop_fitness(i) = (pop_fitness(i)^exponential);
        pop_mem_fitness = 0;%Resets the population member fitness 
    end
    max_fit = max(pop_fitness);%Normalizes the population fitness scores
    pop_fitness = pop_fitness / max_fit;
    fitness_scores = pop_fitness;
end

  