function [fitness_scores] = imagecalculatefitness(target,population,tolerance,exponential,population_size)
    [a b c] = size(target);
    pop_fitness = zeros(population_size,1);
    for i = 1:population_size
        pop_mem = population(:,:,:,i);
        pop_mem_fitness = sum(sum(sum((abs(target - pop_mem))<= tolerance)));
        pop_mem_fitness = pop_mem_fitness / numel(pop_mem);
        pop_fitness(i)= pop_mem_fitness;
        pop_fitness(i) = (pop_fitness(i)^exponential);
        pop_mem_fitness = 0;
    end
    max_fit = max(pop_fitness);
    pop_fitness = pop_fitness / max_fit;
    fitness_scores = pop_fitness;
end

  