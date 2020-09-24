% This function generates n random strings made up of upper case letters,
% lower case letters, and spaces that are the same length as the target
% string
function population = buildpopulation(target_phrase,population_size)
%tic
gen = [32 65:90 97:122];% Creates a vector with values corresponding to spaces, upper case, and lower case letters in the ascii table
population = cell(population_size,1);% Preallocates an empty 200X1 cell array
pop_mem = [];% Preallocates an empty population member vector
 for i = 1:population_size% Makes a counter from 1 to 200
    a = randi(length(gen),1,length(target_phrase)); 
    pop_mem = gen(a);
    pop_mem = char(pop_mem);% Converts all the values in the population member vector to their corresponding characters in the ascii table
    population{i} = pop_mem;% Assigns the current population member to the ith row of the population cell array
    pop_mem = [];% Makes the population member vector empty so that the first population member vector isn't contatenated with the first and again with every iteration
 end
 population = population;
 %toc
end

        