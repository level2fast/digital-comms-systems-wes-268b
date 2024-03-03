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
ofdm_sig_rx_header = (train_seq_length + cp_length*3 + 32);
%% 2.2.a Find the start time of the packet using the cross-correlation
% In OFDM, a training sequence is a predefined set of symbols or patterns 
% that are known both to the transmitter and the receiver. This sequence is
% used for various purposes, including synchronization, channel estimation,
% and equalization
% first add zeros to training sequence to align it with
% 0 bpsk_training_seq(1:26) zeros(1:12) bpsk_training_seq(27:52) 0
bpsk_training_seq_zero_filled = zeros(train_seq_length/2,1);
bpsk_training_seq_zero_filled(2:27) = bpsk_training_seq(1:26);
bpsk_training_seq_zero_filled(39:64) = bpsk_training_seq(27:end);

% modulate the bpsk training seq to an OFDM symbol
bpsk_training_seq_modulated = ifft(bpsk_training_seq_zero_filled);
[rcv_pkt_start, lags] = c_corr(ofdm_sig_rx,bpsk_training_seq_modulated);
x = 0:length(rcv_pkt_start)-1;

[ymax,idx] = max(abs(rcv_pkt_start));
figure(1);
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
%start = idx - (trainer_seq_length/2);
train_seq_start = idx - (train_seq_length/2);
ofdm_sig = ofdm_sig_rx(train_seq_start:end);

train_seq_end = (idx + (train_seq_length/2)-1);
train_seq_samp= train_seq_length/2;

% average over every sample contained with the training sequence to
% determine the frequency offset
train_seq1_start = train_seq_start;
train_seq1_end   = (train_seq_start + fft_size) - 1;
train_seq1       = ofdm_sig_rx(train_seq1_start:train_seq1_end);

train_seq2_start = train_seq1_end + 1;
train_seq2_end   = train_seq2_start + (fft_size-1);
train_seq2       = ofdm_sig_rx(train_seq2_start:train_seq2_end);

z1 = conj(train_seq1) .*  train_seq2;
z1 = mean(z1);
disp(z1)

% determine the frequency offset 
freq_offset = angle(z1) / (2*pi*train_seq_samp);
%% 2.2.a.ii Once you have determined f0, correct your entire packet by applying a frequency offset of −f0 
% (i.e. de-rotating your packet).
%freq_offset_vector = atan2(img(z1)/real(z1))/(2*pi*train_seq_samp);
sig_len = 0:length(ofdm_sig)-1;
ofdm_sig_rx_shifted = exp(-1j*2*pi*freq_offset*sig_len) .* ofdm_sig;
%% 2.2.a.iii After you have removed the frequency offset, estimate the frequency domain channel coefficients H
train_seq = ofdm_sig_rx_shifted(1:train_seq_length);
temp1 = train_seq(1:64);
temp2 = train_seq(65:128);
h_hat1    = fft(temp1,fft_size)   .* (bpsk_training_seq_zero_filled');
h_hat2    = fft(temp2,fft_size) .* (bpsk_training_seq_zero_filled');
h_hat_mat = [h_hat1; h_hat2];
h_hat_avg = (h_hat1 + h_hat2)/2;
h_freqs = linspace(-0.5,0.5,numel(h_hat_avg));
figure(2);
plot(h_freqs,fftshift(abs(h_hat_avg)))
hold on
freqs2 = linspace(-0.5,0.5,numel(ofdm_sig_rx_shifted));
plot(freqs2,abs(fftshift(fft(ofdm_sig_rx_shifted))))
hold off
xlabel('Samples');
ylabel('Magnitude');
title('2.2.a.iii');
%% 2.2.a.iv For the first OFDM block find and remove the cyclic prefix (CP) and then demodulate the OFDM block using the FFT.
header = (train_seq_length + 32 + 16) + 1;
last = header + fft_size;
ofdm_sig_cp_rmvd  = ofdm_sig_rx_shifted(header:last);
ofdm_sig_fft = fft(ofdm_sig_cp_rmvd,fft_size);
figure(3)
scatterplot(ofdm_sig_fft)
title('2.2.a.iv');
%% 2.2.a.v For the first OFDM block, apply the frequency domain equalizer G.
channel_equalization_factor = 1./h_hat_avg;
ofdm_sig_equalized = ofdm_sig_fft .* channel_equalization_factor;
figure(4)
scatterplot(ofdm_sig_equalized)
title('2.2.a.v');

% figure(7)
% figure(8)
%% 2.2.a.vi Apply the equalizer to the entire packet.
process_data = 0;
ofdm_sig_header_rmvd = ofdm_sig_rx_shifted(header:end);
start = header;
ofdm_block_size = 64;
ofdm_sig_final = zeros(fft_size * ofdm_blocks,1)';

% channel_equalization_factor = 1./h_hat_avg;
% ofdm_sig_equalized = ofdm_sig_fft .* channel_equalization_factor;
count = 0;
final_end = 0;
for idxBlks = 1:ofdm_blocks
    startIdx = (idxBlks*ofdm_block_size) +1;
    endIdx =  (startIdx + ofdm_block_size);
    sum = startIdx + endIdx;
    diff = endIdx - startIdx;
    extracted_ofdm_symbol = ofdm_sig_header_rmvd(startIdx:endIdx-1);
    if(idxBlks == 1)
        ofdm_sig_final(idxBlks:fft_size) = fft(extracted_ofdm_symbol,fft_size);
    else
        final_start = (fft_size + final_end) + 1;
        final_end  = (final_start + fft_size)-1;
        fft_out = fft(extracted_ofdm_symbol,fft_size);
        vec = final_start:final_end;
        ofdm_sig_final(vec) = fft_out .* channel_equalization_factor;
    end

end
ofdm_sig_equalized = ofdm_sig_final;
%figure(4)
scatterplot(ofdm_sig_final)
title('2.2.a.vi');

