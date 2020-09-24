clear all;close all;clc

load('asen3300mod.mat')
% sound(signal)
t = linspace(0,(1/fs)*length(signal),length(signal));
figure(1)
plot(t(1:1000),signal(1:1000))

fftsig = fft(signal);
% fftsig = fftshift(fftsig);
f = (fs/2)*linspace(-1,1,length(signal));

figure;
plot(f, abs(fftsig));
title('magnitude FFT of Signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');

fmdemod = demod(signal,fc,fs,'fm');
% sound(fmdemod,fs);

figure
plot(t,fmdemod)

fftsig = fft(fmdemod);
% fftsig = fftshift(fftsig);
f = (fs/2)*linspace(-1,1,length(signal));

figure;
plot(f, abs(fftsig));
title('magnitude FFT of FM Demodulated Signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');

amdemod = demod(signal,fc,fs,'am');
% sound(amdemod,fs)

figure
plot(t,amdemod)

fftsig = fft(amdemod);
% fftsig = fftshift(fftsig);
f = (fs/2)*linspace(-1,1,length(signal));

figure;
plot(f, abs(fftsig));
title('magnitude FFT of AM Demodulated Signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');

load('asen3300mod_noisy.mat')
t = linspace(0,(1/fs)*length(signalnoisy),length(signalnoisy));
plot(t,signalnoisy)
%% analysis 5e
ammod = modulate(fmdemod,fc,fs,'am');

fftsig = fft(ammod);
% fftsig = fftshift(fftsig);
f = (fs/2)*linspace(-1,1,length(signal));

figure;
plot(f, abs(fftsig));
title('magnitude FFT of AM Modulated Signal');
xlabel('Frequency (Hz)');
ylabel('magnitude');
