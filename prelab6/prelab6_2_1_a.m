clear
% Parameters
num_symbols = 128;       % Number of QPSK symbols
samples_per_symbol = 32; % Sampling rate
symbol_rate = 1;         % Symbol rate
fc = 2e6;                % Carrier frequency
fs = samples_per_symbol * symbol_rate; % Sampling frequency

% Generate random binary data for QPSK modulation
data = randi([0 1], 1, 2 * num_symbols);

% Map binary data to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2 * data(1:2:end) - 1 + 1j * (2 * data(2:2:end) - 1));

% Upsample symbols to desired sampling rate
upsampled_symbols = upsample(qpsk_symbols, samples_per_symbol);

% Generate rectangular baseband pulses
rect_pulse = ones(1, samples_per_symbol);

% Pulse shaping using convolution
pulse_shaped_waveform = conv(upsampled_symbols, rect_pulse);

% Modulate the pulse-shaped waveform onto the subcarrier
t = (0:length(pulse_shaped_waveform)-1) / fs; % Time vector
subcarrier = real(pulse_shaped_waveform .* exp(1j * 2 * pi * fc * t));

% Plotting
figure(1);
% Time series plot
subplot(2, 1, 1);
plot(t, subcarrier);
title('QPSK Modulated Subcarrier - Time Series');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Power spectrum plot
subplot(2, 1, 2);
NFFT = 2^nextpow2(length(subcarrier)); % Next power of 2 from length of y
subcarrier_fft = fft(subcarrier,NFFT)/length(subcarrier);
f = fs/2 * linspace(0,1,NFFT/2+1);
subcarrier_fft = subcarrier_fft(1:NFFT/2+1);
psdx = abs(subcarrier_fft).^2;
plot(f,pow2db(psdx));
title('QPSK Modulated Subcarrier - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
grid on;