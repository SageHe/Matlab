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
plot_scope_spectrum_analyzer(2,tvec,ysquare,1e-3, [0 20e3], [-150 0])
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
% computing the autocorrelation
for n = 0:1022
    for i = 0:1022
        sum1(i+1) = CA(i+1)*CA((i+1)+(n));
    end
    Auto(n+1) = 1023/1023 * sum(sum1);
end
%plot the autocorrelation function for maximal length code
figure; stairs(0:1022,Auto); xlim([-10 1030])
title('Autocorrelation Function for Maximal Length Code')
xlabel('Shift/Delay (Chips)');ylabel('Correlation')
y = [CA CA CA CA CA];
Tchip = 1/chiprate;
ind = floor(tvec/Tchip)+1;
ysamp = y(ind);
plot_scope_spectrum_analyzer(3,tvec,ysamp,1e-4, [0 6*chiprate], [-100 0])