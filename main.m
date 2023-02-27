clear
clc

% Read audio files
[xI, ~] = audioread('xI.wav');
[xQ, fs] = audioread('xQ.wav');

% Generate output signal
x = sender(xI, xQ);

% Send signal over channel
y = TSKS10channel(x);

% Recieve signal and interpret 
[zI, zQ, A, tau] = receiver(y);

% Check SignalNoiseRatio, must be greater than 25 dB
SNRzI = 20*log10(norm(xI)/norm(zI-xI));
SNRzQ = 20*log10(norm(xQ)/norm(zQ-xQ));