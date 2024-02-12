% Parameters
num_symbols = 100;      % Number of OFDM symbols
num_carriers = 64;      % Number of subcarriers
cp_length = 16;         % Length of cyclic prefix

% Generate random data for each symbol
data_symbols = randi([0, 1], num_symbols, num_carriers);

% IFFT to convert data to time domain
time_domain_symbols = ifft(data_symbols, [], 2);

% Add cyclic prefix
time_domain_symbols_cp = [time_domain_symbols(:, end - cp_length + 1:end), time_domain_symbols];

% Serialize symbols
ofdm_signal = reshape(time_domain_symbols_cp.', 1, []);

% Plot the generated OFDM signal
figure;
subplot(2, 1, 1);
plot(real(ofdm_signal));
title('Real Part of OFDM Signal');
xlabel('Sample');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(imag(ofdm_signal));
title('Imaginary Part of OFDM Signal');
xlabel('Sample');
ylabel('Amplitude');

% Display some information
disp('OFDM Signal Generation Complete');
disp(['Number of OFDM Symbols: ', num2str(num_symbols)]);
disp(['Number of Subcarriers: ', num2str(num_carriers)]);
disp(['Cyclic Prefix Length: ', num2str(cp_length)]);
