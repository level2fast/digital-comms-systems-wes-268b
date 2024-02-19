%% 2.c use a parallelizer to transform the serial data from a QPSK symbol mapper fed by a random bitstream
% Generate random binary data for QPSK modulation
clear
% Parameters
num_symbols_qpsk = 32;                  % Number of QPSK symbols
samples_per_symbol = 32;               % Sampling rate
rolloff = 0;                           % Rectangular pulse shaping (no roll-off)
symbol_rate = 1;                       % Symbol rate
fc = 2e3;                              % Carrier frequency
fs = samples_per_symbol * symbol_rate; % Sampling frequency

% Generate random binary data for QPSK modulation
data = randi([0 1], 1, 2*num_symbols_qpsk);

% Map binary data to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2 * data(1:2:end) - 1 + 1j * (2 * data(2:2:end) - 1));

% Upsample symbols to desired sampling rate
%upsampled_symbols = upsample(qpsk_symbols, samples_per_symbol);

% Generate rectangular baseband pulses
rect_pulse = ones(1, samples_per_symbol);

% Pulse shaping using convolution
pulse_shaped_waveform = conv(qpsk_symbols, rect_pulse);

% % Modulate the pulse-shaped waveform onto the subcarrier
% t = (0:length(pulse_shaped_waveform)-1) / fs; % Time vector
% subcarrier = pulse_shaped_waveform .* exp(1j * 2 * pi * fc * t);
% 
% % Add a 2 new subcarriers
% N = samples_per_symbol;
% T = symbol_rate;
% % Adjust plot layout
% subplot(2, 1, 1);
% f1_delta  = (1/N*T);
% [~, sub_carrier2] = add_subcarrier(subcarrier,fc + f1_delta,fs,num_symbols_qpsk,N);
% 
% N = samples_per_symbol;
% T = symbol_rate;
% f2_delta  = (1/N*T);
% [ofdm_signal, sub_carrier3] = add_subcarrier(sub_carrier2,fc - f2_delta,fs,num_symbols_qpsk,N);

% % create random bits for input to QPSK symbol mapper
% bitstream = randi([0 1], 1, 2*num_symbols);

% Map binary data to QPSK symbols
%qpsk_symbols = qpsk_symbol_mapper(bit_vector=bitstream,input="bitstream");
s_k = ofdm_parallelizer(symbol_values=ofdm_signal,subcarrier_num=4); % TODO: This isn't working because of size of ofdm_signal
temp  = ifft(s_k);
s_t = ofdm_serializer(symbol_subcarrier_mat=temp);
%% c.i Plot the magnitude of power spectrum of the total transmitted signal s(t) with all subcarriers enabled
% Power spectrum plot
figure(2)
subcarrier = s_t;
subplot(2, 1, 1);
NFFT = 2^nextpow2(length(subcarrier)); % Next power of 2 from length of y
Y = fft(subcarrier,NFFT)/length(subcarrier);
f = fs/2 * linspace(0,1,NFFT/2+1);
ydft = Y(1:NFFT/2+1);
psdx = abs(ydft);
plot(f,pow2db(psdx));
title('QPSK Modulated 4 Subcarriers - Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power(db)');
grid on;
return
%% c.ii zero out subcarriers, Plot the power spectrum of s(t) with only three subcarriers enabled.
% subcarrier = s_t;
% subplot(2, 1, 1);
% NFFT = 2^nextpow2(length(subcarrier)); % Next power of 2 from length of y
% Y = fft(subcarrier,NFFT)/length(subcarrier);
% f = fs/2 * linspace(0,1,NFFT/2+1);
% ydft = Y(1:NFFT/2+1);
% psdx = abs(ydft).^2;
% plot(f,pow2db(psdx));
% title('QPSK Modulated Subcarrier - Power Spectrum');
% xlabel('Frequency (Hz)');
% ylabel('Power(db)');
% grid on;
%% Plot the real part of the corresponding time series for each case. 

% What happens if you try and view all of the subcarriers  constellations 
% overlaid on top of one another? What does this mean if you try and generate 
% an eye-pattern to find the optimal sampling time?






