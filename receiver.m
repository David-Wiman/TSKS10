function [zI, zQ, A, tau] = receiver(y)
% Takes one column vector, x, and returns two signals, zI and zQ, together
% with the delay and scaling of the channel

%% Init
f1 = 75000;
f2 = 95000;
fc = (f1 + f2)/2; % 85 kHz
fs1 = 20000;
fs2 = 400000;
Ts = 1/fs2;
B = 5000;
f_nyq = fs2/2; %found in matlab docs

%% Band-pass filter, degree taken from lab-pm
band_pass_filter = fir1(200, [(f1)/f_nyq (f2)/f_nyq], 'bandpass');
y_filtered = filter(band_pass_filter, 1, y);

% Compute and compensate delay, should equal half the filter order = 100
% https://se.mathworks.com/help/signal/ug/compensate-for-the-delay-introduced-by-an-fir-filter.html
y_filtered = y_filtered(101:end);
y_filtered(end+1:end+100) = zeros(100,1);

%% Demodulate
time_vector = 0:Ts:length(y_filtered)*Ts - Ts;
zI_unfiltered = 2*y_filtered.*cos(2*pi*fc*time_vector)';
zQ_unfiltered = 2*y_filtered.*sin(2*pi*fc*time_vector)';

% Low-pass filter, degree taken from lab-pm, 0.1 given from lab-assisten
low_pass_filter = fir1(200, 0.1, 'low');
zI = filter(low_pass_filter, 1, zI_unfiltered);
zQ = -filter(low_pass_filter, 1, zQ_unfiltered);

% Compute and compensate filter delay
zI = zI(101:end);
zQ = zQ(101:end);
zI(end+1:end+100) = zeros(100,1);
zQ(end+1:end+100) = zeros(100,1);

%% Compute signal delay tau
time_vector = 0:Ts:2-Ts; 
pulseform = chirp(time_vector, 0, 2, 500).*0.2; % Two second chirp
[c_zI, zI_lags] = xcorr(zI, pulseform);
[c_zQ, zQ_lags] = xcorr(zQ, pulseform);
[zI_value, zI_index] =  max(abs(c_zI));
[zQ_value, zQ_index] =  max(abs(c_zQ));
if zI_value >= zQ_value
    tau_index = zI_lags(zI_index);
else
    tau_index = zQ_lags(zQ_index);
end
tau = tau_index*Ts; % 185.00 microseconds

% Compensate signal delay
y_filtered = y_filtered(tau_index+1:end);

%% Demodulate again with signal delay compensation
time_vector = 0:Ts:length(y_filtered)*Ts - Ts;
zI_unfiltered = 2*y_filtered.*cos(2*pi*fc*time_vector)';
zQ_unfiltered = 2*y_filtered.*sin(2*pi*fc*time_vector)';

% Low-pass filter, degree taken from lab-pm, 0.1 given from lab-assisten
low_pass_filter = fir1(200, 0.1, 'low');
zI = filter(low_pass_filter, 1, zI_unfiltered);
zQ = -filter(low_pass_filter, 1, zQ_unfiltered);

% Compensate filter delay
zI = zI(101:end);
zQ = zQ(101:end);
zI(end+1:end+100) = zeros(100,1);
zQ(end+1:end+100) = zeros(100,1);

%% Compute signal scaling
[c_zQ, ~] = xcorr(zQ(1:length(pulseform)), pulseform);
[c_pulse, ~] = xcorr(pulseform);
[~, zQ_index] =  max(abs(c_zQ));
[~, pulse_index] =  max(abs(c_pulse));
nominator = c_zQ(zQ_index);
denominator = c_pulse(pulse_index);
A = nominator/denominator;

% Compensate signal scaling, remove chirp/zeros
zI = zI(length(pulseform)+1:end)/A;
zQ = zQ(length(pulseform)+1:end)/A;

%% Downsample
zI = downsample(zI, 20);
zQ = downsample(zQ, 20);

% Resize to make matrix dimensions match
zI = zI(1:100000);
zQ = zQ(1:100000);


samples = 1:length(zQ);
t = samples/fs1;
plot(t, zQ);
ylabel('Amplitud');
xlabel('Tid [s]');

end