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
active_subc = [1,26,38,63];
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

% modulate the bpsk training seq to an OFDM symbol
bpsk_training_seq_modulated = ifft(bpsk_training_seq);

% N1 is the 1st sample of the 1st training sequence
% N2 is the 1st sample of the second training sequence
% Llong is the number of samples in the each training sequence
tr1_idx = (cp_length * 2) + 1;
tr2_idx = tr1_idx + (train_seq_length/2);
rcv_pkt_start = xcorr(ofdm_sig_rx,bpsk_training_seq_modulated);
disp(max(real(rcv_pkt_start)))
figure;
plot(abs(rcv_pkt_start))
hold on

%% 2.2.a.i Determine the frequency offset f0
% average over every sample contained within each training sequence

%% 2.2.a.ii Once you have determined f0, correct your entire packet by applying a frequency offset of −f0 
% (i.e. de-rotating your packet).

%% 2.2.a.iii After you have removed the frequency offset, estimate the frequency domain channel coefficients H

%% 2.2.a.iv For the first OFDM block find and remove the cyclic prefix (CP) and then demodulate the OFDM block using the FFT.

%% 2.2.a.v For the first OFDM block, apply the frequency domain equalizer G.

%% 2.2.a.vi Apply the equalizer to the entire packet.

