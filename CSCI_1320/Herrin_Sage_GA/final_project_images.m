target = im2double(x36x36monalisa);
population_size = input('How many members are in the population?\n');
population = imagebuildpopulation(target,population_size);
exponential = input('What is the exponential rate?\n');
tolerance = input('What is the tolerance?\n');
%population_fitness = imagecalculatefitness(target,population,tolerance,exponential,population_size);
mating_factor = input('What is the mating factor?\n');
mutationrate = input('What is the mutation rate?\n');
range = input('Whats your desired range of mutation?\n');
randomrate = input('whats the randomrate?\n');
%mating_pool = buildmatingpoolimages(population_fitness,mating_factor);
child = zeros(size(target));
generation = 0;
while generation < 1500
population_fitness = imagecalculatefitness(target,population,tolerance,exponential,population_size);
mating_pool = buildmatingpoolimages(population_fitness,mating_factor);
    for j = 1:population_size
        %population_fitness = calculatefitness(target_phrase,population,population_size);
        %mating_pool = buildmatingpool(population_fitness);
        a = mating_pool(randi(length(mating_pool)));
        par1 = population(:,:,:,a);
        b = mating_pool(randi(length(mating_pool)));
        while b == a
            b = mating_pool(randi(length(mating_pool)));
        end
        par2 = population(:,:,:,b);
        child = imagebreed(par1,par2);
        child = imagemutation(child,mutationrate,range,randomrate);
        population(:,:,:,j) = child;
        a = [];
        b = [];
    end
     generation = generation + 1;
end
child = uint8(child);
imshow(child);
imwrite(child,'Child_image,jpg')

