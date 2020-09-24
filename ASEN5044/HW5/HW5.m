%% Question 2, 2.16
clear;close all;clc

nums = rand(50,1); %Generate 50 random numbers between 0 and 1 from the uniform dist.
figure(1)
histogram(nums,10)
title('Plot of 50 Random Indepented Numbers Between 0 and 1')
xlabel('Bins')
ylabel('Number of Occurences')

mean(nums) % Calculate mean of generated random numbers

std(nums) %Calculate standard deviation of random numbers 

nums500 = rand(500,1);
figure(2)
histogram(nums500,10)
title('Plot of 500 Random Indepented Numbers Between 0 and 1')
xlabel('Bins')
ylabel('Number of Occurences')

mean(nums500)

std(nums500)

nums5000 = rand(5000,1);
figure(3)
histogram(nums5000,10)
title('Plot of 5000 Random Indepented Numbers Between 0 and 1')
xlabel('Bins')
ylabel('Number of Occurences')

mean(nums5000)

std(nums5000)
%% Problem 2.16
clear;close all;clc
for i = 1:10000
    num = -.5 + rand(30,1);
    A(i) = (num(1) + num(2))/2;
    B(i) = (sum(num(1:4)))/4;
    C(i) = (sum(num(1:30)))/30;
end
figure
subplot(1,3,1)
histogram(A,50)
title('(X_1+X_2)/2)')
xlabel('Mean')
subplot(1,3,2)
histogram(B,50)
title('(X_1+X_2+X_3+X_4)/4')
xlabel('Mean')
subplot(1,3,3)
histogram(C,50)
title('(X_1+X_2+...+X_{30})/30')
xlabel('Mean')

