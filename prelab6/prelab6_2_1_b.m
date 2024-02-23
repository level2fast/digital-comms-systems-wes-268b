%% b) Repeat part (a) adding two additional modulated subcarriers spaced ±fd = 1/T =
%% 1/(NTs) apart from the original subcarrier where T is the OFDM block symbol duration.
clear
clf
% Parameters
num_symbols = 128;                     % Number of QPSK symbols
samples_per_symbol = 32;               % Sampling rate
rolloff = 0;                           % Rectangular pulse shaping (no roll-off)
symbol_rate = 1;                       % Symbol rate
fc = 2e3;                              % Carrier frequency
fs = samples_per_symbol * symbol_rate; % Sampling frequency

% Generate random binary data for QPSK modulation
data = randi([0 1], 1, 2*num_symbols);

% Map binary data to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2 * data(1:2:end) - 1 + 1j * (2 * data(2:2:end) - 1));
qpsk_symbols = repelem(qpsk_symbols,samples_per_symbol);

% Modulate the pulse-shaped waveform onto the subcarrier
t = (0:length(qpsk_symbols)-1) / fs; % Time vector
sub_carrier = qpsk_symbols .* exp(1j * 2 * pi * fc * t);

% % Add a 2 new subcarriers
N = samples_per_symbol;
T = symbol_rate;
% Adjust plot layout
Ts = t(end);
f1_delta  = (1/15*N*Ts);
[~, sub_carrier2] = add_subcarrier(sub_carrier,fc + f1_delta,fs,num_symbols,N);

f2_delta  = (1/15*N*Ts);
[ofdm_signal, sub_carrier3] = add_subcarrier(sub_carrier2,fc - f2_delta,fs,num_symbols,N);

figure(2)
subplot(3, 1, 1);
NFFT             = 2^nextpow2(length(sub_carrier)); % Next power of 2 from length of signal
ofdm_fft1_1      = fft(sub_carrier,NFFT)/length(sub_carrier);
ofdm_shifted1_1  = (ofdm_fft1_1);
fshift1_1        = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd1_1      = abs(ofdm_shifted1_1).^2/NFFT; 
plot(fshift1_1, pow2db(ofdm_psd1_1));
title('Single QPSK Modulated Subcarrier(1) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
legend
grid on;

% Power spectrum plot
subplot(3, 1, 2);
NFFT            = 2^nextpow2(length(sub_carrier2)); % Next power of 2 from length of signal
signal1_2       = sub_carrier2;
ofdm_fft1_2     = fft(signal1_2,NFFT)/length(signal1_2);
ofdm_shifted1_2 = fftshift(ofdm_fft1_2);
fshift1_2       = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd1_2     = abs(ofdm_shifted1_2).^2/NFFT; 
hold on
plot(fshift1_1, pow2db(ofdm_psd1_1));
plot(fshift1_2, pow2db(ofdm_psd1_2));
hold off
title('Single QPSK Modulated Subcarrier(1-2) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
legend
grid on;

subplot(3, 1, 3);
NFFT             = 2^nextpow2(length(sub_carrier3)); % Next power of 2 from length of signal
ofdm_fft1_3      = fft(sub_carrier3,NFFT)/length(sub_carrier3);
ofdm_shifted1_3  = fftshift(ofdm_fft1_3);
fshift1_3        = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd1_3      = abs(ofdm_shifted1_3).^2/NFFT; 
hold on
plot(fshift1_1, pow2db(ofdm_psd1_1));
plot(fshift1_2, pow2db(ofdm_psd1_2));
plot(fshift1_3, pow2db(ofdm_psd1_3));
hold off
title('Single QPSK Modulated Subcarrier(1-3) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
legend
grid on;