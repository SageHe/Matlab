function [out] = anagram1(word)
wordan = '';%Make wordan an empty string
vec = 1:length(word);%Creates a vector "vec" of values from one to the number of letters in the word input in to the function
for a = 1:length(word)%Creates a counter from one to the lenght of the input word
    disp(vec)%Displays the created vector
    b = vec(randi(numel(vec)));%creates a variable b that is vec of some random number from inside vec
    disp(b);%Displays this variable b
    c = find(vec == b);%Finds where vec is equivalent to the variable b found earlier 
    vec(c) = [];%Deletes the value represented by c from the original vector
    wordan(a) = word(b);%Makes wordan of the value of a equal to the value of word of b
    out = (wordan);%Makes the output of the function equal to the string of wordan
end
end