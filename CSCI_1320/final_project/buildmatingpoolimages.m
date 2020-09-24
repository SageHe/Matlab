%Image mating pool
function out = buildmatingpool(pop_fitness,mating_factor)
%tic
mating_pool = []; %Initializes mating pool to be empty vector
%mating_factor = input('What is your mating factor?\n'); %Saves the rate of the mating factor as a variable
    for i = 1:length(pop_fitness) %Counter goes from one to the length of the population
        num_tickets = (pop_fitness(i) * mating_factor); %The number of tickets for that population member is the rounded number of the fitness of that member multiplied by the mating factor
        mating_pool = [mating_pool (i*ones(1,round(num_tickets)))]; % Adds the number of indexes of the ith pop member to the mating pool
    end
out = mating_pool; %Output the mating pool
%toc
end        
