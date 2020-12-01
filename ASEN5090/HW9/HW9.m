close all;clc
%% Part 1 -- Visibility Prediction
% read in given almanac data
C =  2.99792458e8; % speed of light,m/s
data = load('HW9data.mat');
[gps_alm,gps_alm_cell] = read_GPSyuma('YUMA310.ALM',2);
userpos_ecef = lla2ecef([40.01043021 -105.24388889 1614.976]); %negative 105 due to W not E

weeknum = 2130; %week 82 without modulo 1024
tow = 591660;
prn = [1:14 16:32];
for i = prn
    [health(i),pos(i,:)] = broadcast_eph2pos(gps_alm,[weeknum tow],i);
    [az(i),el(i),range(i)] = compute_azelrange(userpos_ecef,pos(i,:));
end
%% Part 2 -- Carrier Wipeoff/Code Correlation
Sr = data.gpsdata;

tvec = linspace(0,1e-3,5001);
prn31 = CA_Gen(31,1024);
% prn31 = flip(prn31);
prn31 = (prn31==0)*1 + (prn31==1)*(-1);
chiprate = 1.023e6;
Tchip = 1/chiprate;
ind = floor(tvec/Tchip)+1;
sampprn = prn31(ind);

IF  = 1.25e6; %intermediate freq.
fd = 0; %Doppler fequency 

% theta = 2*pi*(IF+fd)*tvec;
delay = 9;
n = numel(tvec);
S = comp_corr(Sr,delay,tvec,sampprn,IF,fd);
%% Part 3 -- Delay/Doppler Grid Search
delayax = 0:1:4999; %change this to 5001 and see if it still works
doppax = -10e3:500:10e3;
S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn,IF,doppax(j)));
    end
end    
%% Section for plotting mesh plot    
[X,Y] = meshgrid(delayax,doppax);

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 4
% found 32, 1, 22, 10, repeat part 3 above for prns 32, 1, and 22
% Starting with PRN 1 for part 4
vis = find(el > 10);
prn1 = CA_Gen(vis(1),1024);
prn1 = (prn1==0)*1 + (prn1==1)*(-1);
sampprn1 = prn1(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn1,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 4 -- PRN 22
prn22 = CA_Gen(vis(6),1024);
prn22 = (prn22==0)*1 + (prn22==1)*(-1);
sampprn22 = prn22(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn22,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 4 -- PRN 32
prn32 = CA_Gen(vis(9),1024);
prn32 = (prn32==0)*1 + (prn32==1)*(-1);
sampprn32 = prn32(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn32,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 4 -- PRN 10
prn10 = CA_Gen(vis(3),1024);
prn10 = (prn10==0)*1 + (prn10==1)*(-1);
sampprn10 = prn10(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn10,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 5 -- Increase integration time and repeat
tvec = linspace(0,2e-3,10001);
ind = floor(tvec/Tchip)+1;
prn1 = CA_Gen(1,2048);
prn1 = (prn1==0)*1 + (prn1==1)*(-1);
sampprn1 = prn1(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn1,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 5 -- PRN 10
prn10 = CA_Gen(10,2048);
prn10 = (prn10==0)*1 + (prn10==1)*(-1);
sampprn10 = prn10(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn10,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor
%% Part 5 -- PRN 22
prn22 = CA_Gen(22,2048);
prn22 = (prn22==0)*1 + (prn22==1)*(-1);
sampprn22 = prn22(ind);

S = zeros(numel(delayax),numel(doppax));
for i = 1:numel(delayax)
    for j = 1:numel(doppax)
        S(i,j) = abs(comp_corr(Sr,delayax(i),tvec,sampprn22,IF,doppax(j)));
    end
end    

normS = S./max(S(:));
normS = normS';


[peakdop,peakdel] = find(normS == 1);

fprintf('Doppler shift: %0.2f Hz \n',doppax(peakdop))
fprintf('Delay: %0.0f samples \n',delayax(peakdel))

figure
surf(X,Y,normS)
xlabel('Phase Shift (samples)','FontSize',15)
ylabel('Doppler Shift (Hz)','FontSize',15)

dop = normS(peakdop,:);
del = normS(:,peakdel);
dt = (tvec(2) - tvec(1));

figure
plot(delayax.*dt,dop)
xlabel('\tau (s)')
ylabel('S')
grid on
grid minor

figure
plot(doppax,del)
xlabel('Doppler (Hz)')
ylabel('S')
grid on
grid minor