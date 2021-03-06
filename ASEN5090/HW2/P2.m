clear all;close all;clc
%% Question 2.1, part a -- Create entire 1023-chip C/A code for PRN 19 as a
%vector. Plot first 16 and last 16 chips of code. Express in hex, check
%that first 16 chips of PRN 19 are E6D6 in hex.
%Start with G1 and G2 registers as all 1's
G1 = ones(1,10);
G2 = ones(1,10);
%PRN 19 uses phase selection 3 & 6
for i = 1:1023
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(3),G2(6));
    CA_PRN_19(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
figure
subplot(1,2,1)
plot(CA_PRN_19(1:16),'*')
title('First 16 chips of PRN 19')
xlabel('Chip Number')
ylabel('Chip Value')
subplot(1,2,2)
plot([1008:1023],CA_PRN_19(1008:end),'*')
title('Last 16 chips of PRN 16')
xlabel('Chip Number')
ylabel('Chip Value')
%Part b -- do agin for epoch 1024-2046 and compare 1023-element vector or one
%from part a (should be exact same?)
for i = 1:1023
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(3),G2(6));
    CA2(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
%Part c -- repeat part a for PRN 25 -> use phase selection of 5 and 7
G1 = ones(1,10);
G2 = ones(1,10);
for i = 1:1023
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(5),G2(7));
    CA_PRN_25(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
figure
subplot(1,2,1)
plot(CA_PRN_25(1:16),'*')
title('First 16 chips of PRN 25')
xlabel('Chip Number')
ylabel('Chip Value')
subplot(1,2,2)
plot([1008:1023],CA_PRN_25(1008:end),'*')
title('Last 16 chips of PRN 25')
xlabel('Chip Number')
ylabel('Chip Value')
%Part d -- repeat part a for PRN 5 -> use phase selection of 1 and 5
G1 = ones(1,10);
G2 = ones(1,10);
for i = 1:1023
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(1),G2(9));
    CA_PRN_5(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
figure
subplot(1,2,1)
plot(CA_PRN_5(1:16),'*')
title('First 16 chips of PRN 5')
subplot(1,2,2)
plot([1008:1023],CA_PRN_5(1008:end),'*')
title('Last 16 chips of PRN 5')
%% Question 2.2
%Part a
%Begin by converting CA codes for given PRNs from previous problem from
%{0,1} -> {1,-1} via 0 -> 1 and 1 -> -1
G1 = ones(1,10);
G2 = ones(1,10);
%PRN 19 uses phase selection 3 & 6
for i = 1:4096
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(3),G2(6));
    CA_PRN_19(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
CA_PRN_25_C = (CA_PRN_25 == 0)*1 + (CA_PRN_25 == 1)*-1;
CA_PRN_5_C = (CA_PRN_5 == 0)*1 + (CA_PRN_5 == 1)*1;
%Generate auto-correlation function of PRN 19
figure
hold on
for n = 0:2046
    for i = 1:1023
        this(i) = CA_PRN_19_C(i)*CA_PRN_19_C(i+n);
    end
    R_19(n+1) = (1/1023)*sum(this);
%     plot(n,R_19,'o','--b')
end
plot([1:numel(R_19)],R_19)
%% part b -- Create 1023-chip PRN sequence as PRN 19 delayed by 200 chips and
%plot normalized cylic cross-correlation of delayed sequence with original
%PRN 19 sequence
clear CA_PRN_19 CA_PRN_19_c
G1 = ones(1,10);
G2 = ones(1,10);
%PRN 19 uses phase selection 3 & 6
for i = 1:4092
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
    PS = xor(G2(3),G2(6));
    CA_PRN_19(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
%shift PRN 19 by 200 chips to make alternate code
clear this
n = 200;
CA_Replica = CA_PRN_19_C;
for i = 1:n
    CA_Replica = CA_Replica([end (1:end - 1)]);
end
%Cross-correlation of original and shifted PRN 19
for j = 1:1023
    for k = 1:1023
        this(k) = CA_PRN_19_C(k)*CA_Replica(k+j);
    end
    R_lk(j) = (1/1023)*sum(this);
end
figure
plot([1:numel(R_lk)], R_lk)
clear CA_PRN_19 CA_PRN_19_C CA_PRN_25 CA_PRN_25_C
%% Part c -- Plot cyclic cross-correlation of PRN 19 and PRN 25
CA_PRN_19 = CA_Gen(19,1023);
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
CA_PRN_25 = CA_Gen(25,1023);
CA_PRN_25_C = (CA_PRN_25 == 0)*1 + (CA_PRN_25 == 1)*-1;
n = length(CA_PRN_19_C);
x = reshape(CA_PRN_19_C,1,n);
yshifted = reshape(CA_PRN_25_C,n,1);
lag = [0:n]; Rxy = NaN(size(lag));
for i = 1:length(lag)
    Rxy(i) = x*yshifted*(1/1023);
    yshifted = [yshifted(end); yshifted(1:end-1)];
end
% for n = 1:1023
%     for i = 1:1023
%         this(i) = CA_PRN_19_C(i)*CA_PRN_25_C(i+n);
%     end
%     R_lk(n) = (1/1023)*sum(this);
% end
figure
plot(lag,Rxy,'.',lag,Rxy,'-'),grid,xlim([0 n]);
%% Part d -- Plot cyclic cross-correlation of PRN 19 and PRN 5
CA_PRN_19 = CA_Gen(19,1023);
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
CA_PRN_5 = CA_Gen(5,1023);
CA_PRN_5_C = (CA_PRN_25 == 0)*1 + (CA_PRN_25 == 1)*-1;
n = length(CA_PRN_19_C);
x = reshape(CA_PRN_19_C,1,n);
yshifted = reshape(CA_PRN_5_C,n,1);
lag = [0:n]; Rxy = NaN(size(lag));
for i = 1:length(lag)
    Rxy(i) = x*yshifted*(1/1023);
    yshifted = [yshifted(end); yshifted(1:end-1)];
end
figure
plot(lag,Rxy,'.',lag,Rxy,'-'),grid,xlim([0 n]);
%% Part e -- Sum given 3 PRNs together and correlate with PRN 19
CA_PRN_19 = CA_Gen(19,2046);
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
x1 = CA_Gen(19,2046);
x1 = (x1 == 0)*1 + (x1 == 1)*-1;
x1 = [x1(end-349:end) x1(1:end-350)];
x2 = CA_Gen(25,2046);
x2 = (x2 == 0)*1 + (x2 == 1)*-1;
x2 = [x2(end-904:end) x2(1:end-905)];
x3 = CA_Gen(5,2046);
x3 = (x3 == 0)*1 + (x3 == 1)*-1;
x3 = [x3(end-74:end) x3(1:end-75)];
x = x1+x2+x3;
n = length(x);
y = reshape(CA_PRN_19_C,1,n); 
yshifted = reshape(y,n,1);
lag = [-length(x)/2:length(x)/2]; Rxy = NaN(size(lag));
for i = 1:length(lag)
    Rxy(i) = x*yshifted*(1/1023);
    yshifted = [yshifted(end); yshifted(1:end-1)];
end
figure
plot(lag,Rxy,'.',lag,Rxy,'-'),grid,xlim([-n/2 n/2]);
% title('
xlabel('Shift Number')
ylabel('Normalized Cyclic Correlation')
%% Part f -- Plot x1,x2,x3, and noise all with same vertical axis scale
x1 = CA_Gen(19,1023);
x1 = (x1 == 0)*1 + (x1 == 1)*-1;
x1 = [x1(end-349:end) x1(1:end-350)];
x2 = CA_Gen(25,1023);
x2 = (x2 == 0)*1 + (x2 == 1)*-1;
x2 = [x2(end-904:end) x2(1:end-905)];
x3 = CA_Gen(5,1023);
x3 = (x3 == 0)*1 + (x3 == 1)*-1;
x3 = [x3(end-74:end) x3(1:end-75)];
noise = 4*randn(1,1023); %Generate noise based on given mean and standard deviation
figure
subplot(2,2,1)
plot([0:length(x1)-1],x1),grid,ylim([-15 15])
title('PRN 19 350-chip Shift')
xlabel('Shift Number')
ylabel('PRN Value')
subplot(2,2,2)
plot([0:length(x2)-1],x2),grid,ylim([-15 15])
title('PRN 25 905-chip Shift')
xlabel('Shift Number')
ylabel('PRN Value')
subplot(2,2,3)
plot([0:length(x3)-1],x3),grid,ylim([-15 15])
title('PRN 5 75-chip Shift')
xlabel('Shift Number')
ylabel('PRN Value')
subplot(2,2,4)
plot([0:length(noise)-1],noise),grid,ylim([-15 15])
xlabel('Shift Number')
ylabel('PRN Value')
%% Part g -- Sum x1,x2,x3, and noise and correlate with original PRN 19
x = x1+x2+x3+noise;
CA_PRN_19 = CA_Gen(19,1023);
CA_PRN_19_C = (CA_PRN_19 == 0)*1 + (CA_PRN_19 == 1)*-1;
n = length(x);
y = reshape(CA_PRN_19_C,1,n); 
yshifted = reshape(y,n,1);
lag = [0:n]; Rxy = NaN(size(lag));
for i = 1:length(lag)
    Rxy(i) = x*yshifted*(1/1023);
    yshifted = [yshifted(end); yshifted(1:end-1)];
end
figure
plot(lag,Rxy,'.',lag,Rxy,'-'),grid,xlim([0 n]);
title('Correlation VS Shifts')
xlabel('Shift Number')
ylabel('Normalized Cyclic Correlation')




