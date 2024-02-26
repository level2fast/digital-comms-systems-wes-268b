%% 1.b Now add cyclic prefix (CP) to the data that is longer than the delay in the channel model.
% Summary: In OFDM, the cyclic prefix (CP) is a guard interval inserted at 
% the beginning of each symbol to mitigate the effects of multipath 
% propagation and inter-symbol interference (ISI) in wireless communication
% channels.
% 
% Benefits:
% The CP isolates different OFDM blocks from each other when the 
% wireless channel contains multiple paths, i.e. is frequency-selective.
% The CP turns the linear convolution with the channel into a circular 
% convolution. Only with a circular convolution, we can use the single-tap 
% equalization OFDM is so famous for.
%
% Negative effects in OFDM:
% Requires more data when sending an OFDM signal.
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
num_symbols_qpsk    = 512;                    % Number of QPSK symbols
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

% add a cyclic prefix to tx signal
% -------------------------------------------------------------------------
% the CP of an OFDM symbol is obtained by prepending a copy of the last NCP
% samples from the end of the OFDM signal to its beginning. This way we 
% obtain a circular signal structure, i.e. the first NCP and last NCP 
% samples are equal in each OFDM symbol.
cp_length = 22;
ofdm_num_subc = num_subcarriers + cp_length;
ofdm_num_symb = size(temp_normalized,2);
ofdm_norm_cp = zeros(ofdm_num_subc,ofdm_num_symb);
for symbTxIdx = 1:ofdm_num_symb
    ofdm_norm_cp(:,symbTxIdx) = add_cyclic_prefix(temp_normalized(:,symbTxIdx),cp_length); 
end 

% convert symbols from parallel to serial
sig_serial = ofdm_serializer(symbol_subcarrier_mat=ofdm_norm_cp);

% create transmit signal
sig_tx = sig_serial;

% need to create a delayed version of the received signal. This simulates the
% effect that the channel impulse response has on the transmitted OFDM signal.
%
% Convolution can be used to create a delayed and shifted version of a
% signal. Convolve the transmitted signal with a discrete time impulse
% response.

% create the channel impulse response
channel_imp_rsp = create_ch_imp_rsp(taps=32,time_delay=16);

% convolve transmit signal with channel impulse response
sig_rx = conv(sig_tx,channel_imp_rsp);

%% Repeat part prelab7.2.a
% figure(1)
% shifted_sig     = sig_rx(1:length(sig_tx));
% %shifted_sig     = circshift(sig_tx,16);
% subcarrier_all  = shifted_sig;
% NFFT            = 2^nextpow2(length(subcarrier_all)); % Next power of 2 from length of signal
% signal1_1       = subcarrier_all;
% ofdm_fft1_1     = fft(signal1_1,NFFT);
% fshift1_1       = (-NFFT/2:NFFT/2-1);
% ofdm_psd1_1     = abs(ofdm_fft1_1).^2; 
% plot(fshift1_1, fftshift(pow2db(ofdm_psd1_1)));
% title('2.1.a.i OFDM RX Signal w/ Channel Impulse Response - PSD');
% xlabel('Frequency (Hz)');
% ylabel('Power(dB)');
% %ylim([-80 80])
% grid on;

%% Plot the constellation of several subcarriers 
shifted_sig     = sig_rx(1:length(sig_tx));
%shifted_sig     = circshift(sig_tx,16);

% working solution
% cp_length = 22;
% cp_segments = temp_normalized((end-cp_length+1):end,:);
% extended_x = [cp_segments; temp_normalized];
% serial_x = reshape(extended_x,[],1);
% shifted_sig = circshift(serial_x,16);
% shifted_sig = serial_x + shifted_sig;

shifted_sig_parallel     = reshape(shifted_sig, ofdm_num_subc, []);
ofdm_symbols_mat = zeros(num_subcarriers,ofdm_num_symb);
for symbRxIdx = 1:ofdm_num_symb
    ofdm_symbols_mat(:,symbRxIdx) = remove_cyclic_prefix(shifted_sig_parallel(:,symbRxIdx),cp_length); 
end 
recv_offset_sig = fft(ofdm_symbols_mat,num_subcarriers);

% working solution
% faded_x = reshape(shifted_sig,size(extended_x));
% fft_in = faded_x(((end-num_subcarriers+1):end),:);
% recv_offset_sig = fft(fft_in,num_subcarriers);

figure(2)
scatter(real(recv_offset_sig(1,:)),imag(recv_offset_sig(1,:)))
title("Constellation of 1st subcarrier");
ylabel("Imaginary")
xlabel("Real")

figure(3)
scatter(real(recv_offset_sig(2,:)), imag(recv_offset_sig(2,:)))
title("Constellation of 2nd Subcarrier");
ylabel("Imaginary")
xlabel("Real")

figure(4)
scatter(real(recv_offset_sig(3,:)), imag(recv_offset_sig(3,:)))
title("Constellation of 3rd Subcarrier");
ylabel("Imaginary")
xlabel("Real")