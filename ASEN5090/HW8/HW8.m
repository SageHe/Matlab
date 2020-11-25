clear alll;close all;clc
%% Sampling and Frequency Resolution
% a.) create time vec. for sampling frequency of 50 MHz, starting at t=0
% with duration of 10 msec.
f_samp = 50e6;
tvec = [0:(1/f_samp):10e-3];
tvec(end) = [];
%b.) nyquist frequency is at max half the sampling frequency -> 25 MHz
f_nyq = 2*f_samp;
f_res = f_samp/numel(tvec);
% duration needed to produce frequency res. of 5 Hz is 200 msec

%% Oscilloscope and Spectrum Analyzer Plots
%a.) create 200 KHz sine wave w/ amp. of 1V. Plot signal as seen on oscope
%and spectrum analyzer
y = sin(2*pi*200e3*tvec);

plot_scope_spectrum_analyzer(1,tvec,y,1e-5, [0 400e3], [-400 0])
%b.) describe key features of both plots, how spectrum analyzer result
%compares with what would be expected based off of fourier series rep.

%% Square Wave
%a.) Create 10 KHz square wave with amp. of 1V. Plot signal as it would be
%seen on an oscope and on a spectrum analyzer
ysquare = square(2*pi*10e3*tvec);
plot_scope_spectrum_analyzer(2,tvec,ysquare,1e-3, [0 2e5], [-150 0])
%b.) describe key features of both plots, how spectrum analyzer result
%compares to what would be expected based on fourier series rep. of square wave

%% PRN Codes - Maximal Length
%a.) create maximal length prn code using GPS C/A code G1 register with
%chipping rate of 
chiprate = 1.023e6;
G1 = ones(1,10);
XGi = G1(10);
for i = 1:2045
    % update G1 
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    % the maximal length PRN Code
    XGi = [G1(10) XGi];
end
% the maximal length C/A Code
CA = (XGi==0)*1 + (XGi==1)*(-1);
CA = flip(CA);
% computing the autocorrelation
for n = 0:1022
    for i = 0:1022
        sum1(i+1) = CA(i+1)*CA((i+1)+(n));
    end
    Auto(n+1) = 1023/1023 * sum(sum1);
end
%plot the autocorrelation function for maximal length code
% figure; stairs(0:1022,Auto); xlim([-10 1030])
% title('Autocorrelation Function for Maximal Length Code')
% xlabel('Shift/Delay (Chips)');ylabel('Correlation')
y = [CA CA CA CA CA];
% y = flip(y);
Tchip = 1/chiprate;
ind = floor(tvec/Tchip)+1;
ysamp = y(ind);
figure
stairs(tvec(485:495),ysamp(485:495))
plot_scope_spectrum_analyzer(4,tvec,ysamp,1e-4, [0 6*chiprate], [-100 0])
figure
stairs(tvec(485:495),ysamp(485:495),'LineWidth',2); 
ylim([-1.2 1.2])
xlabel('Time (s)')
ylabel('Chip Value')
title('Chip Value VS Time for 485-495 indices') 

% xlim([-10 1030])
%% PRN Codes - Gold Code
%a.) create GPS C/A code for PRN 5 with chipping rate of 1.023 MHz
PRN5 = CA_Gen(5,2046);
PRN5 = (PRN5==0)*1 + (PRN5==1)*-1;
for n = 0:1022
    for i = 0:1022
        sum2(i+1) = PRN5(i+1)*PRN5((i+1)+(n));
    end
    Auto2(n+1) = 1023/1023 * sum(sum2);
end
%plot the autocorrelation function for maximal length code
% figure; stairs(0:1022,Auto2); xlim([-10 1030])
% title('Autocorrelation Function for PRN5 Code')
% xlabel('Shift/Delay (Chips)');ylabel('Correlation')
PRN5 = flip(PRN5);
PRN5 = [PRN5 PRN5 PRN5 PRN5 PRN5];
prn5samp = PRN5(ind);
figure
stairs(tvec(485:495),prn5samp(485:495),'LineWidth',2); 
ylim([-1.2 1.2])
xlabel('Time (s)')
ylabel('Chip Value')
title('Chip Value VS Time for 485-495 indices') 
% xlim([-10 1030])
plot_scope_spectrum_analyzer(8,tvec,prn5samp,1e-4, [0 6*chiprate], [-100 0])

%% Direct Spread Spectrum Modulation
car_sig = cos(2*pi*tvec*5*chiprate);
mod_sig = car_sig.*ysamp;
plot_scope_spectrum_analyzer(9,tvec,mod_sig,2e-6, [0 10*chiprate], [-100 0])

%% BOC codes
%a.) implement BOC(1,1) code using G1 for PRN code and same carrier, code
%chipping rate, and time
BOC11 = square(pi*tvec*2*chiprate);
bocmod_sig = BOC11.*mod_sig;
%b.) plot signal as it would be seen on an oscilloscope and on a spectrum
%analyzer
plot_scope_spectrum_analyzer(10,tvec,bocmod_sig,2e-6, [0 10*chiprate], [-100 0])
%c.) describe key features of both plots. Discuss width of main lobe and
%side loves and identify the code repeat frequency  based on the spectrum
%analyzer results 

%% Demodulation/Carrier Recovery
%a.) Add white noise with standard deviation of 1V to the modulated signal
%created in 6a.
noise = normrnd(0,1,[1,numel(ysamp)]);
noise_sig = mod_sig + noise;
plot_scope_spectrum_analyzer(11,tvec,noise_sig,2e-6, [0 10*chiprate], [-100 0])
noisy_sig_rec = noise_sig.*ysamp;
plot_scope_spectrum_analyzer(12,tvec,noisy_sig_rec,2e-6, [0 10*chiprate], [-100 0])
