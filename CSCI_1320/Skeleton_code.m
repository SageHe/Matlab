%Task 1.1
    %Generate a random string of characters the same length as the target
    %phrase by using a for loop with i from 1 to the length(target) creating
    %random values corresponding with the ascii table values for
    %uppercase letters, lowercase letters, and spaces, changing those
    %randomely generated values to characters using char(), then putting
    %them in a vector. After every random string is generated, put that random string
    %in a structure so that each random string can be referenced with just an index


%Task 1.2  
    %Compare each member of the population to the target phrase by using a
    %for loop with a counter from one to the number of characters in the
    %target phrase and a logical statement that evaluates if pop.member(i)
    %is equivalent to target(i) and then adds the returned logical value to
    %a vector. Sum the vector and divide that sum by either the number of
    %numbers in the pop.member vector or the number of characters in the
    %target phrase in order to calculate the fitness of each population
    %member. Multiply this number by 100 to express as a percentage to
    %"normalize" the fitness score.
    
    % add up all fitness scores and divide by the number of test populations
    % this will calculate the average fitness of all test populations as a whole

%Task 1.3
	%use a for loop to compare the fitness score of every test population and find the
	%maximum score. 
%Divide all scores by this value in order to “normalize”
	%In order to assign an appropriate number of “mating raffle tickets”, multiply
	%each normalized score by 10 and round (using round() function) 
	%Note: Each “ticket” should simply be an index for which population is to be 
	% used in mating, not the entire string

%Task 1.4
	%First method
	%Generate random number between 1 and length of target string
	%”Child” string will be created using two for loops; first from one to random 
	%number, with each “child” char the same as parent 1
	%Second for loop from rand number+1 to length, each “child” char same as 
	%parent 2
	%Second method
	%Create vector with same length as target, random integers 0 to 1
	%Create opposite vector, first by making zeros vector with same length
	%use find(~ ) function to assign value of 1 to indexes w/ value of 0 in first vector 
	%child gets first parent characters from “1” locations of first vector, second parent
	%characters from “1” locations of second vector

%Task 1.5
	%Define mutation rate as decimal
	%Multiply rate by 100
	%Create (length of target) random integers from 1 to modified rate
	%If random number is equal to 1, assign character a new random character

%Task 1.6 
	%Receive input of target phrase, input into function from Task 1.1
	%output should be array made up of each target and each test population
	%Use Task 1.2 function to calculate fitness scores for each test population
	%Call Task 1.3 function to build mating pool
	%Select “parent” populations using raffle system
	%Breed using Task 1.4 function
	%Ask for input of mutation rate
	%Use this in Task 1.5 function in order to apply random mutations to populations








%Task 2.1
	%Recieve input (target) image, and separate into three channels (RGB)
	%Create test populations with same (3D) dimensions as test image
	%Assign random numbers equal to pixel intensity (0-255)
	
%Task 2.2
%Method 1
%Using nested for loops, compare each pixel of each test population to that of
	%the target image
	%For each test pixel within 10 of the target value, add 1 to fitness score of that 
%test population 
%Method 2
%Use mean filter function from Assignment 7
%Method 3
%Use diff function to look at how fast change is happening from pixel to pixel (in
% 2 dimensions)
%If a pixel is changing very fast relative to its neighbor, this area most likely has
%low fitness
%Experiment to combine these three methods in order to find an efficient and %effective method

%Task 2.3
    % Randomly change the brightness of a pixel instead of completely
    % randomly mutating it by saving a mutation range variable and instead
    % of mutating each pixel to a new number in the original range (0-255),
    % instead change its value by adding a random number in the stated
    % mutation range. For example, to mutate it by a set random certain
    % amount, if you desired mutation range is 30, create a range of
    % numbers from (-30,30) and add a random number from this range to the
    % current value of the pixel. Do this by using a nested for loop that
    % goes through the image and adds a new random value in the stated
    % mutation range to every pixel
    
    %In order to make it a 1/4 chance that the pixel will be completely
    %randomly mutated instead of in the desired mutation range, create a
    %vector of values from one to four. Instate a if else statement that
    %randomly picks a value from the vector of values from one to four, and
    %if that selected value is equal to one of the four values in that
    %vector, mutate the pixel in the full range (0-255), else, use the
    %origianl mutation range and method stated above.

%Task 2.4
    %The image was originally manipulated in color in task 2.1, but if this
    %is not done and it is originally used in black and white, at the end
    %of every mutation of the complete image, combine each layer of R G % B
    %in order to make it a color image instead of just a single-layer black
    %and white image
