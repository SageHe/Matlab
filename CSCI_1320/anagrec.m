function [scrambled] = anagrec(word)
if length(word) == 1
    anagram = word;
else
a = randi(numel(word));%Picks a random number in the range from one to the number of elements in the word and assigns it to the variable a
b = word(a);%Assigns the ath letter of the word to the variable b
word(a) = [];%Deletes that corresponding letter from the original word
anagram = [b];%Stores the chosen random letter from the word in the variable anagram
scrambled = [anagram, anagrec(word)];%Concatenates the anagram along with calling the function again
end
end


