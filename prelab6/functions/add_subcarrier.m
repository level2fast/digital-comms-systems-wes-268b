function [output_signal,new_subcarrier] = add_subcarrier(signal,f_delta,fs,num_symbols,samples_per_symbol)
%ADD_SUBCARRIER Summary of this function goes here
%   Detailed explanation goes here

% Generate random binary data for QPSK modulation
data = randi([0 1], 1, 2*num_symbols);

% Map binary data to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2 * data(1:2:end) - 1 + 1j * (2 * data(2:2:end) - 1));
qpsk_symbols = repelem(qpsk_symbols,samples_per_symbol);

% Add subcarrier by performing a frequency shift
t = (0:length(qpsk_symbols)-1) / fs; % Time vector
new_subcarrier = exp(1j * 2 * pi * f_delta * t) .* (qpsk_symbols);
output_signal = signal .* new_subcarrier;
end

