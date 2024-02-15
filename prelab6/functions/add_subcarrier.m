function [output_signal,new_subcarrier] = add_subcarrier(signal,f_delta,fs,num_symbols,samples_per_symbol)
%ADD_SUBCARRIER Summary of this function goes here
%   Detailed explanation goes here

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
new_subcarrier = pulse_shaped_waveform .* exp(1j * 2 * pi * f_delta * t);
output_signal = signal .* new_subcarrier;
end

