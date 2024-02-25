%% 1.a Use the 32 subcarrier data from the previous simulation problem in Prelab 6
% Summary: Delay spread refers to the spreading of a transmitted signal over 
% time due to multipath propagation in the channel. As a result, different 
% copies of the transmitted signal arrive at the receiver with different 
% delays, causing distortion of the received signal. The delay spread can 
% be characterized by the spread of the impulse response of the channel.
%
% Negative effects in OFDM:
% Delay spread can cause intersymbol interference (ISI) and degrade 
% system performance.
%
% How to model:
% One common model used to describe delay spread in OFDM systems is the 
% tapped delay line model. In this model, the channel impulse response is 
% represented as a sum of delayed and attenuated copies of the transmitted 
% signal, referred to as taps. Each tap corresponds to a distinct path in
% the multipath channel.
clear
clf
% Simulate the channel given by the discrete impulse response h[n] = δ[n] − δ[n − 16].
% create the received signal
num_symbols_qpsk    = 128;                    % Number of QPSK symbols
num_bits_per_symbol = 2;
total_qpsk_symbols  = num_bits_per_symbol*num_symbols_qpsk;
num_subcarriers     = 32;
num_samp_per_symbol = 32;
num_samples = total_qpsk_symbols * num_samp_per_symbol;

% Generate random binary data for QPSK modulation
data = randi([0 3],[total_qpsk_symbols * num_subcarriers 1] );

% Map binary data to QPSK symbols;
qpsk_symbols = exp(1j*data*pi/2);
sub_carrier = qpsk_symbols;

% convert symbols from serial to parallel
sig_parallel = ofdm_parallelizer(symbol_values=sub_carrier);

% perform ifft to on parallelized data 
temp_normalized  = ifft(sig_parallel,num_subcarriers) * sqrt(num_subcarriers);

% convert symbols from parallel to serial
sig_serial = ofdm_serializer(symbol_subcarrier_mat=temp_normalized);

% create transmit signal
sig_tx = sig_serial;

% need to create a delayed version of the received signal. This simulates the
% effect that the channel impulse response has on the transmitted OFDM signal.

% Convolution can be used to create a delayed and shifted version of a
% signal. Convolve the transmitted signal with a discrete time impulse
% response.

% create the channel impulse response
channel_imp_rsp = create_ch_imp_rsp(taps=16,time_delay=16);

% convert received signal to parallel

% iterate over each subcarrier index
for subcIdx = 1:32
    % perform convolution with CIR on each subcarrier
    sig_rx = conv(sig_tx,channel_imp_rsp);
end

%% 1.a.i Plot the power density spectrum for the received signal given this channel model
figure(1)
shifted_sig = circshift(sig_tx,16);
subcarrier_all = shifted_sig + sig_tx;
subplot(2, 1, 1);
NFFT            = 2^nextpow2(length(subcarrier_all)); % Next power of 2 from length of signal
signal1_1       = subcarrier_all;
ofdm_fft1_1     = fft(signal1_1,NFFT);
fshift1_1       = (-NFFT/2:NFFT/2-1);
ofdm_psd1_1     = abs(ofdm_fft1_1).^2; 
plot(fshift1_1, fftshift(pow2db(ofdm_psd1_1)));
title('2.1.a.i OFDM RX Signal w/ Channel Impulse Response - PSD');
xlabel('Frequency (Hz)');
ylabel('Power(dB)');
grid on;

%% 1.a.ii Plot the constellation of several subcarriers 
%figure(6)
shifted_sig = circshift(subcarrier_all,16);
shifted_sig = reshape(shifted_sig, num_subcarriers, []);

recv_offset_sig = fft(shifted_sig);

figure('name',['Matlab Simulation 1.a.ii']);
scatter(real(recv_offset_sig(1,:)),imag(recv_offset_sig(1,:)))
title("Constellation of the Lowest Frequency Subcarrier w/ a 2 ..." + ...
    "Sample Offset");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1.a.ii']);
scatter(real(recv_offset_sig(num_subcarriers/4,:)), imag(recv_offset_sig(num_subcarriers/4,:)))
title("Constellation of the Middle Subcarrier w/ a 2 sample offset");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1.a.ii']);
scatter(real(recv_offset_sig(num_subcarriers/2,:)), imag(recv_offset_sig(num_subcarriers/2,:))         )
title("Constellation of the Highest Frequency Subcarrier w/ a 2 ..." + ...
    "sample offset");
ylabel("Imaginary")
xlabel("Real")


%% 1.b Now add cyclic prefix(CP) to the data that is longer than the delay in the channel model