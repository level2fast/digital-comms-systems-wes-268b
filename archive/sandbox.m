clear
clf
% % Parameters
% num_symbols = 8; % Number of QPSK symbols
% samples_per_symbol = 32; % Sampling rate
% rolloff = 0; % Rectangular pulse shaping (no roll-off)
% symbol_rate = 1; % Symbol rate
% fc = 2e6; % Carrier frequency
% fs = samples_per_symbol * symbol_rate; % Sampling frequency
% 
% % Generate random binary data for QPSK modulation
% data = randi([0 1], 1, 2*num_symbols);
% 
% % Map binary data to QPSK symbols
% qpsk_symbols = 1/sqrt(2) * (2 * data(1:2:end) - 1 + 1i * (2 * data(2:2:end) - 1));
% 
% % Upsample symbols to desired sampling rate
% upsampled_symbols = upsample(qpsk_symbols, samples_per_symbol);
% 
% % Generate rectangular baseband pulses
% rect_pulse = ones(1, samples_per_symbol);
% 
% % Pulse shaping using convolution
% pulse_shaped_waveform = conv(upsampled_symbols, rect_pulse);
% 
% % Modulate the pulse-shaped waveform onto the subcarrier
% t = (0:length(pulse_shaped_waveform)-1) / fs; % Time vector
% subcarrier = real(pulse_shaped_waveform .* exp(1i * 2 * pi * fc * t));
% 
% % Plotting
% figure;
% 

% % Time series plot
% subplot(2, 1, 1);
% plot(t, subcarrier);
% title('QPSK Modulated Subcarrier - Time Series');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % Power spectrum plot
% subplot(2, 1, 2);
% NFFT = 2^nextpow2(length(subcarrier)); % Next power of 2 from length of y
% Y = fft(subcarrier,NFFT)/length(subcarrier);
% f = fs/2*linspace(0,1,NFFT/2+1);
% plot(f,2*abs(Y(1:NFFT/2+1)));
% title('QPSK Modulated Subcarrier - Power Spectrum');
% xlabel('Frequency (Hz)');
% ylabel('Power');
% grid on;
% 
% % Adjust plot layout
% subplot(2, 1, 1);

% % Define parameters
% num_symbols = 32;  % Number of QPSK symbols
% rolloff_factor = 0.5; % Roll-off factor for pulse shaping
% 
% % Generate random binary data for QPSK modulation
% binary_data = randi([0 1], num_symbols, 2); % 2 bits per symbol
% 
% % Map binary data to QPSK symbols
% qpsk_symbols = 2 * binary_data - 1; % Map 0 to -1, 1 to 1
% 
% % Perform QPSK modulation
% qpsk_modulated = qpsk_symbols(:,1) + 1i * qpsk_symbols(:,2);
% 
% % Plot QPSK signal on a scatter plot
% figure;
% scatter(real(qpsk_modulated), imag(qpsk_modulated), 'filled');
% xlabel('In-phase');
% ylabel('Quadrature');
% title('QPSK Signal Constellation');
% axis([-1.5 1.5 -1.5 1.5]); % Set axis limits for better visualization
% grid on;

% % Define parameters
% N = 1000; % Number of symbols
% fc = 1e6; % Carrier frequency (1 MHz)
% fs = 10e6; % Sampling frequency (10 MHz)
% T = 1/fs; % Sampling period
% t = (0:N-1)*T; % Time vector
% SNR_dB = 20; % Signal-to-Noise Ratio in dB
% 
% % Generate random binary data for QPSK modulation
% data = randi([0 1], 1, N);
% 
% % Modulation: QPSK
% s = 1/sqrt(2) * (2*data - 1) + 1j/sqrt(2) * (2*randi([0 1], 1, N) - 1); % QPSK symbols
% 
% % Upconvert to passband
% x = s .* exp(1j * 2 * pi * fc * t);
% 
% % Add Gaussian noise
% noise_power = 10^(-SNR_dB/10);
% noise = sqrt(noise_power/2) * (randn(1, N) + 1j*randn(1, N));
% received_signal = x + noise;
% 
% % Downconvert to baseband
% received_signal = received_signal .* exp(-1j * 2 * pi * fc * t);
% 
% % Demodulation: QPSK
% demodulated_data = (real(received_signal) > 0) + 1j*(imag(received_signal) > 0);
% 
% % Calculate Bit Error Rate (BER)
% errors = sum(data ~= (real(demodulated_data) > 0));
% BER = errors/N;
% 
% disp(['Bit Error Rate (BER): ', num2str(BER)]);
% QPSK symbol generation using complex exponentials
% 



% % Number of symbols
% num_symbols = 1000;
% 
% % Generate random binary data for QPSK modulation
% data = randi([0 1], 1, 2*num_symbols);
% 
% % Convert binary data to QPSK symbols
% qpsk_symbols = zeros(1, num_symbols);
% 
% for i = 1:num_symbols
%     % Map binary pairs to complex symbols
%     if data(2*i-1) == 0 && data(2*i) == 0
%         qpsk_symbols(i) = 1/sqrt(2) * exp(1i*pi/4); % 00 -> (1 + j)
%     elseif data(2*i-1) == 0 && data(2*i) == 1
%         qpsk_symbols(i) = 1/sqrt(2) * exp(1i*3*pi/4); % 01 -> (-1 + j)
%     elseif data(2*i-1) == 1 && data(2*i) == 0
%         qpsk_symbols(i) = 1/sqrt(2) * exp(1i*7*pi/4); % 10 -> (1 - j)
%     else
%         qpsk_symbols(i) = 1/sqrt(2) * exp(1i*5*pi/4); % 11 -> (-1 - j)
%     end
% end
% 
% % Plot the constellation diagram of QPSK symbols
% scatterplot(qpsk_symbols);
% xlabel('In-phase');
% ylabel('Quadrature');
% title('QPSK Constellation Diagram');
% 
% % Display generated QPSK symbols
% disp('Generated QPSK symbols:');
% disp(qpsk_symbols);
% Define the number of symbols
num_symbols = 32;

% Generate random bits for QPSK modulation
bits = randi([0, 1], 1, 2*num_symbols);

% Convert bits to symbols
symbols = bi2de(reshape(bits, 2, []).', 'left-msb');

% Map symbols to constellation points
constellation = exp(1j * pi/4 * (2*symbols + 1));

% Plot constellation diagram
figure;
scatterplot(constellation);
title('QPSK Constellation Diagram');
xlabel('In-phase (I)');
ylabel('Quadrature (Q)');
grid on;
axis equal;

% time_shift_us = 1/1e6;   % shift time by 10 microseconds
% t0 = t * (time_shift_us); 
% f  = 1e3;
% omega = 2*pi*f;
% shifted_sig = s_t .* exp(1j*t0*omega);