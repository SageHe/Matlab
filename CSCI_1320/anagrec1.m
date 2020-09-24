function [out] = anagrec1(word)
persistent anagram%Causes the variable anagram to exist through every iteration of the function
if length(word) == 1
    anagram = [anagram word];
    out = anagram;%If the length of the word is one letter, this concatenates the anagram with the remaining letter, produces the created anagram as the out, then clears anagram so it does not exist for the next initial run of the function 
    clear anagram
else
a = randi(numel(word));%Picks a random number in the range from one to the number of elements in the word and assigns it to the variable a
b = word(a);%Assigns the ath letter of the word to the variable b
word(a) = [];%Deletes that corresponding letter from the original word
anagram = [anagram b];%Stores the chosen random letter from the word in the variable anagram
out = anagrec1(word);%Concatenates the anagram along with calling the function again
end
end
