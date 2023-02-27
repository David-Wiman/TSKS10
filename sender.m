function x = sender(xI, xQ)
% Takes two column vectors, xI and xQ, and returns an I/Q modulated signal x

%% Init
f1 = 75000;
f2 = 95000;
fc = (f1 + f2)/2; % 85 kHz
fs1 = 20000;
fs2 = 400000;
Ts = 1/fs2;
B = 5000;
f_nyq = fs2/2; %found in matlab docs

%% Upsample, compensate for enery loss, 20*20 kHz = 400 kHz
xI_upsample = upsample(xI, 20)*20;
xQ_upsample = upsample(xQ, 20)*20;

% Low-pass filter, interpolating
% https://se.mathworks.com/help/signal/ref/fir1.html#bulla52-Wn
low_pass_filter = fir1(200, 2*B*0.917/f_nyq, 'low'); % 0.917 from trail and error
xI_upsample_filtered = filter(low_pass_filter, 1, xI_upsample);
xQ_upsample_filtered = filter(low_pass_filter, 1, xQ_upsample);

% Compute and compensate delay, should equal half the filter order = 100
% https://se.mathworks.com/help/signal/ug/compensate-for-the-delay-introduced-by-an-fir-filter.html
xI_upsample_filtered = xI_upsample_filtered(101:end);
xQ_upsample_filtered = xQ_upsample_filtered(101:end);
xI_upsample_filtered(end+1:end+100) = zeros(100,1);
xQ_upsample_filtered(end+1:end+100) = zeros(100,1);

%% Add known pulse form, chirp
time_vector = 0:Ts:2-Ts; % Two second chirp
pulseform = chirp(time_vector, 0, 2, 500).*0.2; % Make sure chirp amplitude matches signal
xI_upsample_filtered_zeros = [zeros(length(pulseform),1); xI_upsample_filtered];
xQ_upsample_filtered_chirp = [pulseform'; xQ_upsample_filtered]; % xQ gets the chirp

%% Generate output signal
time_vector = 0:Ts:(length(xI_upsample_filtered_zeros))*Ts - Ts;
x = xI_upsample_filtered_zeros.*cos(2*pi*fc*time_vector)' - xQ_upsample_filtered_chirp.*sin(2*pi*fc*time_vector)'; % Peak around 0.175

end