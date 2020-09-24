% This function calculates the fitness of the overall population
function out = calculatefitness(target_phrase,population,population_size,exponential_factor)
%tic
pop_fitness = cell(population_size,1);% Creates an empty cell array that is the population size rows long by one column
pop_mem_fitness = [];%Makes the pop_mem_fitness vector empty
for i = 1:population_size%Creates a counter from one to the population size
    pop_mem = population{i};
    pop_mem_fitness = sum(pop_mem == target_phrase);
    pop_mem_fitness = (pop_mem_fitness / length(pop_mem));
    pop_fitness{i,1} = pop_mem_fitness;% Adds the current value of pop_mem_fitness to the ith cell of the cell array
    pop_mem_fitness = [];% Makes the pop_mem_fitness vector empty so the current value of pop_mem_fitness is not added to the current calculated value
end
pop_fit_max = max([pop_fitness{:}]);
for k = 1:length(pop_fitness)
    pop_fitness{k} = ([pop_fitness{k}])^exponential_factor;
   % pop_fitness{k} = (pop_fitness{k} / pop_fit_max);    
end
pop_fit_max = max([pop_fitness{:}]);
for j = 1:length(pop_fitness)
    pop_fitness{j} = (pop_fitness{j} / pop_fit_max);
end
out = pop_fitness;
%toc
end
            
            
        