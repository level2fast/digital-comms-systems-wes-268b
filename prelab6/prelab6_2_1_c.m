%% 2.c use a parallelizer to transform the serial data from a QPSK symbol mapper fed by a random bitstream
% Generate random binary data for QPSK modulation
clear
clf
% Parameters
num_symbols = 512;                    % Number of QPSK symbols
num_subcarriers  = 32;
num_samp_per_symbol = 32;
num_samples = 2*num_symbols * num_samp_per_symbol;

% Generate random binary data for QPSK modulation
data = randi([0 3],[2*num_symbols*num_subcarriers 1] );

%% Map binary data to QPSK symbols;
qpsk_symbols = exp(1j*data*pi/2);
sub_carrier = qpsk_symbols;
%% c.i Plot the magnitude of power spectrum of the total transmitted signal s(t) with all subcarriers enabled
% convert symbols from serial to parallel
s_k = ofdm_parallelizer(symbol_values=sub_carrier);

% perform ifft to on parallelized data 
temp_normalized  = ifft(s_k,num_subcarriers) * sqrt(num_subcarriers);

% convert symbols from parallel to serial
s_t = ofdm_serializer(symbol_subcarrier_mat=temp_normalized);

% Power spectrum plot
% figure(1)
% subcarrier_all = s_t;
% subplot(2, 1, 1);
% NFFT            = 2^nextpow2(length(subcarrier_all)); % Next power of 2 from length of signal
% signal1_1       = subcarrier_all;
% ofdm_fft1_1     = fft(signal1_1,NFFT);
% fshift1_1       = (-NFFT/2:NFFT/2-1);
% ofdm_psd1_1     = abs(ofdm_fft1_1); 
% plot(fshift1_1, fftshift(mag2db(ofdm_psd1_1)));
% title('2.1.c.i OFDM All Subcarriers - Magnitude');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;


%% c.ii Plot the power spectrum of s(t) with only three subcarriers enabled.
% convert data from serial to parallel
s_k2 = ofdm_parallelizer(symbol_values=sub_carrier);

% grab 3 sub carriers
subcarrier1_3_mat = s_k2;
subcarrier1_3_mat(4:end,:) = 0;
subcarrier1_3_mat_mod = ifft(subcarrier1_3_mat) * sqrt(num_subcarriers);
subcarrier1_3_serial = ofdm_serializer(symbol_subcarrier_mat = subcarrier1_3_mat_mod);

% subplot(2, 1, 2);
% freqs = linspace(-0.5,0.5,numel(s_t));
% signal1_1       = subcarrier1_3_serial;
% ofdm_fft1_1     = abs(fft(signal1_1));
% ofdm_fft1_1_db  = mag2db(ofdm_fft1_1);
% ofdm_shifted1_1 = fftshift(ofdm_fft1_1_db);
% plot(freqs, ofdm_shifted1_1);
% title('2.1.c.ii OFDM 3 Subcarriers - Power Specturm');
% xlabel('Frequency (Hz)');
% ylabel('Power(dB)');
% grid on;

%% Plot the real part of the corresponding time series for each case. 
samp_freq_hz = 1e3;
samp_period_s = 1/samp_freq_hz;
total_dur_s = (num_samples-1) * samp_period_s;
t = 0:samp_period_s:total_dur_s; 

figure(2)
subcarrier_all = s_t;
t_all = (0:length(sub_carrier)-1); % Time vector
subplot(2, 1, 1);
plot(t, real(subcarrier_all));
title('2.1.c.iii OFDM All Subcarriers - Real Component');
xlabel('Time');
ylabel('Amplitude');
grid on;

subcarrier1_3 = subcarrier1_3_serial;
t1_3 = (0:length(subcarrier1_3)-1); % Time vector
subplot(2, 1, 2);
plot(t, real(subcarrier1_3));
title('2.1.c.iii OFDM 3 Subcarriers - Real Component');
xlabel('Time');
ylabel('Amplitdue');
grid on;

%% Plot the constellation
figure(3)
subplot(2,1,1)
%figure('name','Matlab Simulation 1ciii Constellation of 32 Subcarriers');
scatter(real(subcarrier_all), imag(subcarrier_all))
title("2.1.c.iii Constellation of 32 Subcarriers");
ylabel("Imaginary")
xlabel("Real")

subplot(2,1,2)
%figure('name','Matlab Simulation 1ciii Constellation of 3 Subcarriers');
scatter(real(subcarrier1_3), ...
    imag(subcarrier1_3))
title("2.1.c.iii Constellation of 3 Subcarriers");
ylabel("Imaginary")
xlabel("Real")

%% Plot an eye diagram
% figure(4)
% eyediagram(subcarrier_all,2*num_samp_per_symbol)

%% d-f
figure(5)
% d. generate noise and add it to serialized signal s_t
mean = 0;
std = 0.05;
in_phase_noise = normrnd(mean, std, 1, num_samples);
quad_noise = normrnd(mean, std, 1, num_samples);
noise = in_phase_noise + 1j*quad_noise;
noisy_sig_time = s_t + noise;

% e. Take forward fft of noisy signal
noisy_sig_freq = ofdm_parallelizer(symbol_values=noisy_sig_time);
noisy_sig_freq = fft(noisy_sig_freq);

% f. Generate multiple OFDM frames at least 1000 and plot the QPSK
% constellation
sub_carrier_lowest_freq = noisy_sig_freq(1,:);

subplot(3,1,1)
scatter(real(sub_carrier_lowest_freq),imag(sub_carrier_lowest_freq))
title("Constellation of the Lowest Frequency Subcarrier");
ylabel("Imaginary")
xlabel("Real")

sub_carrier_highest_freq = noisy_sig_freq(num_subcarriers/2,:);
subplot(3,1,2)
scatter(real(sub_carrier_highest_freq),imag(sub_carrier_highest_freq))
title("Constellation of the highest Frequency Subcarrier");
ylabel("Imaginary")
xlabel("Real")

sub_carrier_halfway_freq = noisy_sig_freq(num_subcarriers/4,:);
subplot(3,1,3)
scatter(real(sub_carrier_halfway_freq),imag(sub_carrier_halfway_freq))
title("Constellation of the Halfway Frequency Subcarrier");
ylabel("Imaginary")
xlabel("Real")

% g
figure(6)
shifted_sig = circshift(s_t,1);
shifted_sig = reshape(shifted_sig, num_subcarriers, []);

recv_offset_sig = fft(shifted_sig);

subplot(3,1,1)
scatter(real(recv_offset_sig(1,:)),imag(recv_offset_sig(1,:)))
title("Constellation of the Lowest Frequency Subcarrier w/ a 2 ..." + ...
    "Sample Offset");
ylabel("Imaginary")
xlabel("Real")

subplot(3,1,2)
scatter(real(recv_offset_sig(num_subcarriers/4,:)), imag(recv_offset_sig(num_subcarriers/4,:)))
title("Constellation of the Middle Subcarrier w/ a 2 sample offset");
ylabel("Imaginary")
xlabel("Real")

subplot(3,1,3)
scatter(real(recv_offset_sig(num_subcarriers/2,:)), imag(recv_offset_sig(num_subcarriers/2,:))         )
title("Constellation of the Highest Frequency Subcarrier w/ a 2 ..." + ...
    "sample offset");
ylabel("Imaginary")
xlabel("Real")

