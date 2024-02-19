%% b) Repeat part (a) adding two additional modulated subcarriers spaced ±fd = 1/T =
%% 1/(NTs) apart from the original subcarrier where T is the OFDM block symbol duration.

% Parameters
num_symbols = 256;                     % Number of QPSK symbols
samples_per_symbol = 32;               % Sampling rate
rolloff = 0;                           % Rectangular pulse shaping (no roll-off)
symbol_rate = 1;                       % Symbol rate
fc = 2e3;                              % Carrier frequency
fs = samples_per_symbol * symbol_rate; % Sampling frequency

% Generate random binary data for QPSK modulation
data = randi([0 1], 1, 2*num_symbols);

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
subcarrier = pulse_shaped_waveform .* exp(1j * 2 * pi * fc * t);

% Add a 2 new subcarriers
N = samples_per_symbol;
T = symbol_rate;
% Adjust plot layout

f1_delta  = (1/N*T);
[~, sub_carrier2] = add_subcarrier(subcarrier,fc + f1_delta,fs,num_symbols,N);

f2_delta  = (1/N*T);
[ofdm_signal, sub_carrier3] = add_subcarrier(sub_carrier2,fc - f2_delta,fs,num_symbols,N);

% Plotting
figure(1);
% Time series plot
subplot(2, 1, 1);
plot(t, ofdm_signal);
title('Single QPSK Modulated Subcarriers(1-3) - Time Series');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Power spectrum plot
subplot(2, 1, 2);
NFFT = 2^nextpow2(length(ofdm_signal)); % Next power of 2 from length of signal
ofdm_fft = fft(ofdm_signal,NFFT)/length(ofdm_signal);
ofdm_shifted = fftshift(ofdm_fft);
fshift = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd = abs(ofdm_shifted).^2/NFFT; 
plot(fshift, pow2db(ofdm_psd));
title('Single QPSK Modulated Subcarrier(1-3) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
legend
grid on;



