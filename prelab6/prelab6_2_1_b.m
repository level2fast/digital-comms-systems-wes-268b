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

% Modulate the pulse-shaped waveform onto the subcarrier
t = (0:length(qpsk_symbols)-1) / fs; % Time vector
sub_carrier = qpsk_symbols .* exp(1j * 2 * pi * fc * t);

% Add a 2 new subcarriers
N = samples_per_symbol;
T = symbol_rate;
% Adjust plot layout
Ts = t(end);
f1_delta  = (1/N*Ts);
[~, sub_carrier2] = add_subcarrier(sub_carrier,fc + f1_delta,fs,num_symbols);

f2_delta  = (1/N*Ts);
[ofdm_signal, sub_carrier3] = add_subcarrier(sub_carrier2,fc - f2_delta,fs,num_symbols);

% Plotting
figure(1);
% Time series plot
subplot(3, 1, 1);
NFFT            = 2^nextpow2(length(sub_carrier)); % Next power of 2 from length of signal
signal1_1       = sub_carrier;
ofdm_fft1_1     = fft(signal1_1,NFFT)/length(signal1_1);
ofdm_shifted1_1 = fftshift(ofdm_fft1_1);
fshift1_1       = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd1_1     = abs(ofdm_shifted1_1).^2/NFFT; 
plot(fshift1_1, pow2db(ofdm_psd1_1));
title('Single QPSK Modulated Subcarrier(1) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
hold on
legend
grid on;

% Power spectrum plot
subplot(3, 1, 2);
NFFT            = 2^nextpow2(length(sub_carrier)); % Next power of 2 from length of signal
signal2_1       = sub_carrier;
ofdm_fft2_1     = fft(signal2_1,NFFT)/length(signal2_1);
ofdm_shifted2_1 = fftshift(ofdm_fft2_1);
fshift2_1       = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd2_1     = abs(ofdm_shifted2_1).^2/NFFT; 
plot(fshift2_1, pow2db(ofdm_psd2_1));
title('Single QPSK Modulated Subcarrier(1-2) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
hold on
legend
grid on;

NFFT             = 2^nextpow2(length(sub_carrier2)); % Next power of 2 from length of signal
ofdm_fft2_2      = fft(sub_carrier2,NFFT)/length(sub_carrier2);
ofdm_shifted2_2  = fftshift(ofdm_fft2_2);
fshift2_2        = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd2_2      = abs(ofdm_shifted2_2).^2/NFFT; 
plot(fshift2_2, pow2db(ofdm_psd2_2));
hold on

% Power spectrum plot
subplot(3, 1, 3);
NFFT          = 2^nextpow2(length(sub_carrier)); % Next power of 2 from length of signal
signal3_1       = sub_carrier;
ofdm_fft3_1     = fft(signal2_1,NFFT)/length(signal2_1);
ofdm_shifted3_1 = fftshift(ofdm_fft3_1);
fshift3_1       = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd3_1     = abs(ofdm_shifted3_1).^2/NFFT; 
plot(fshift3_1, pow2db(ofdm_psd3_1));
title('Single QPSK Modulated Subcarrier(1-3) - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
hold on
legend
grid on;

NFFT          = 2^nextpow2(length(sub_carrier2)); % Next power of 2 from length of signal
ofdm_fft3_2     = fft(sub_carrier2,NFFT)/length(sub_carrier2);
ofdm_shifted3_2 = fftshift(ofdm_fft3_2);
fshift3_2        = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd3_2      = abs(ofdm_shifted3_2).^2/NFFT; 
plot(fshift3_2, pow2db(ofdm_psd3_2));
hold on

NFFT          = 2^nextpow2(length(sub_carrier3)); % Next power of 2 from length of signal
ofdm_fft3_3     = fft(sub_carrier3,NFFT)/length(sub_carrier3);
ofdm_shifted3_3 = fftshift(ofdm_fft3_3);
fshift3_3        = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
ofdm_psd3_3      = abs(ofdm_shifted3_3).^2/NFFT; 
plot(fshift3_3, pow2db(ofdm_psd3_3));
hold on