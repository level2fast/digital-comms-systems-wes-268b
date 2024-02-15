%% Prelab6 Problem 2.2.f
% Define the number of samples per bit
N = 8;
% Generate a random bit sequence
bit_sequence = [1 0 0 1 1]; % Generating a sequence of 10 bits

% Modulate the bit sequence using BPSK
% Map 0 to -1 and 1 to 1
bpsk_signal = 2 * bit_sequence - 1;

% Flatten the array
bpsk_signal = bpsk_signal(:)'; 

% Time vector for plotting
t = 1:length(bpsk_signal);

% Plotting the modulated BPSK waveform
figure(1);
subplot(2, 1, 1);
stairs(t, bpsk_signal);
title('Modulated BPSK Waveform');
xlabel('Time');
ylabel('Amplitude');
axis([1 length(bpsk_signal) -1.5 1.5]);
grid on;

% IFFT BPSK using the F^-1 matrix
[idft_twiddle_factors_mat] = calc_idft_twiddle_factors([1 1 1 1 1 1 1 1]);
temp = zeros(8,5);
temp(2,:) = bpsk_signal;
modulated_bpsk =  idft_twiddle_factors_mat * temp ;
moduladated_bpsk_vec = reshape(modulated_bpsk,1,[]);

subplot(2, 1, 2);
plot(real(moduladated_bpsk_vec));
title('Real Part  Modulated BPSK Waveform');
xlabel('Samples');
ylabel('Amplitude');
grid on;

