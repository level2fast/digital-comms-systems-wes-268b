%% 2 Channel Estimation and Equalization
% Summary: 
% Channel estimation in (OFDM): 
% Ensures high-quality data transmission over
% wireless channels. 
%
% Channel Equalization: 
% Once the channel response is estimated for all subcarriers, this 
% information is used to equalize the received signal, effectively mitigating
% the channel impairments. The equalization adjusts the amplitude and phase
% of each subcarrier to counteract the channel's effects.
%
% Benefits: 
% Enables the receiver to accurately decode the transmitted data despite 
% the impairments introduced by the wireless channel 
clear
clf
ofdm_packet = load("ofdm_pkt.mat");
fft_size = 64;
trainer_seq_size = 64;
active_subc = [1:26;38:63];
cp_length = 16;
train_seq_length = 128;
num_bit_per_symbol = 1;
num_symbols = 0;
ofdm_blocks = 10;
bpsk_training_seq = ofdm_packet.s;
ofdm_sig_rx = ofdm_packet.y;
%% 2.2.a Find the start time of the packet using the cross-correlation
% In OFDM, a training sequence is a predefined set of symbols or patterns 
% that are known both to the transmitter and the receiver. This sequence is
% used for various purposes, including synchronization, channel estimation,
% and equalization

% first add zeros to training sequence to align it with
% 0 bpsk_training_seq(1:26) zeros(1:12) bpsk_training_seq(27:52) 0
bpsk_training_seq_zero_filled = [0];
bpsk_training_seq_zero_filled(end + 26)= 0;
bpsk_training_seq_zero_filled(2:27) = bpsk_training_seq(1:26);
bpsk_training_seq_zero_filled(28:37) = 0;
bpsk_training_seq_zero_filled(38:63) = bpsk_training_seq(27:52);
bpsk_training_seq_zero_filled(64) = 0;

% modulate the bpsk training seq to an OFDM symbol
bpsk_training_seq_modulated = ifft(bpsk_training_seq_zero_filled);
%display(bpsk_training_seq_modulated)

[rcv_pkt_start, lags] = c_corr(ofdm_sig_rx,bpsk_training_seq_modulated);
x = 0:length(rcv_pkt_start)-1;
disp(max(abs(rcv_pkt_start)))
figure(1);
[ymax,idx] = max(abs(rcv_pkt_start));
plot(x,abs(rcv_pkt_start))
hold on
plot(x(idx), ymax, 'ro')
hold off
xlabel('Samples');
ylabel('Magnitude)');
title('2.1.a Cross Correlation of training seq with rx signal- Magnitude');

%% 2.2.a.i Determine the frequency offset f0
% average over every sample contained within each training sequence
% N1 is the 1st sample of the 1st training sequence
% N2 is the 1st sample of the second training sequence
% Llong is the number of samples in the each training sequence
% extract the training sequence from the received sequence
train_seq_start = idx - train_seq_length;
train_seq_end = idx;
train_seq = ofdm_sig_rx(train_seq_start:train_seq_end);
train_seq_samp= train_seq_length;

% average over every sample contained with the training sequence to
% determine the frequency offset
% calculate 
z1=0;
for idxTrainSeq = 1:length(train_seq)
    z_1 = z1+ conj(train_seq(idxTrainSeq)) .* train_seq(idxTrainSeq+1);
end
z1 = (1/train_seq_length) * z1;



% determine the frequency offset 


freq_offset_vector = atan2(img(Z2)/real(Z2))/(2*pi*train_seq_samp);

tr1_idx = 129;
tr2_idx = 129 + trainer_seq_size;

%% 2.2.a.ii Once you have determined f0, correct your entire packet by applying a frequency offset of −f0 
% (i.e. de-rotating your packet).

%% 2.2.a.iii After you have removed the frequency offset, estimate the frequency domain channel coefficients H

%% 2.2.a.iv For the first OFDM block find and remove the cyclic prefix (CP) and then demodulate the OFDM block using the FFT.

%% 2.2.a.v For the first OFDM block, apply the frequency domain equalizer G.

%% 2.2.a.vi Apply the equalizer to the entire packet.

