%% Sampling and Frequency Resolution
% a.) create time vec. for sampling frequency of 50 MHz, starting at t=0
% with duration of 10 msec.
tvec = [0:(1/50e6):10e-3];
%b.) nyquist frequency is at minimum double the sampling frequency -> 100 MHz