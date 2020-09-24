tic
max_fitness = [];
avg_fitness = [];
best_phrase = [];
%clear all;close all;clc
global mating_factor
global mutation_rate
target_phrase = input('Enter your target phrase:','s');
population_size = input('How many members are in the population?\n');
population = buildpopulation(target_phrase,population_size);
generation = 1;
child = zeros(1,length(target_phrase));
%child = [];
%for i = 1:length(target_phrase)
    %child = [child 0];
%end
mating_factor = input('What is your mating factor?\n');
exponential_factor = input('What is your exponentail factor?\n');
mutation_rate = input('What is the mutation rate(percentage)?\n');
while sum(child == target_phrase) ~= length(target_phrase)
%for k = 1:20
    population_fitness = calculatefitness(target_phrase,population,population_size,exponential_factor);
    mating_pool = buildmatingpool(population_fitness,mating_factor);
        for j = 1:length(population)
            %population_fitness = calculatefitness(target_phrase,population,population_size);
            %mating_pool = buildmatingpool(population_fitness);
            a = mating_pool(randi(length(mating_pool)));
            par1 = population{a};
            b = mating_pool(randi(length(mating_pool)));
            while b == a
                b = mating_pool(randi(length(mating_pool)));
            end
            par2 = population{b};
            child = breed(par1,par2);
            child = causeMutation(child);
            population{j} = child;
            a = [];
            b = [];
        end
    %end
 max_fitness = [max_fitness; max([population_fitness{:}])];
 avg_fitness = [avg_fitness; (sum([population_fitness{:}])) / length([population_fitness{:}])];
% best_phrase = [best_phrase; population{find(max(population_fitness))}];
generation = generation + 1;
%Write neccessary stuff to text file
%hold on;
%plot(generation, avg_fitness);
end
gen_diversity = (max_fitness - avg_fitness);
write_out = [max_fitness avg_fitness gen_diversity best_phrase];
dlmwrite('write_out.txt',write_out);
fprintf('%dgenerations',generation)      
generation = 1:(generation - 1);
hold on;
plot(generation, avg_fitness);
plot(generation, max_fitness);
toc
%fprintf('%dgenerations',generation)      
    
