%Sage Herrin, 106071909, Assignment 8, 3/25/17, section 104
[num, txt, raw] = xlsread('Section9_data');%Reads the the text, numerical, and raw data from the Excel file into matlab
add = 0;%Makes add equal to zero
total = [];%Preallocates the total matrix to be empty
for i = 3:103%Makes a counter from 3 to 103
    lump_sum = [];%Preallocates lump_sum matrix to be empty
    for j = 7:23%Makes a counter from 7 to 23
        add = (raw{i,j} == raw{2,j});%Assigns add to be equal to the logical statement 1 or 0 of raw{i,j} equivalent to raw{2,j}
        lump_sum = [lump_sum add];%Contatenates lump_sum with itself and the value of add for every iteration in the loop
    end
    score = sum(lump_sum);%Assigns score to be equal to the sum of the values in the lump_sum vector 
    score = (score * 5);%Multiplies score by 5
    score = (score + raw{i,5});%Adds the value of the value in the ith row of the 5th column of every row, or the short answer score
    total = [total score];%Concatenates the value of total and the value of score in to a vector
    clear lump_sum%Clears the lump_sum variable for the next iteration of the ith counter loop
end
histogram(total,10)%Makes a histogram of the calculated scores with a scale of ten on the x axis and the number of students who recieved that score on the y axis
        
        
    