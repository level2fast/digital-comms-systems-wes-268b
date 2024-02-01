%% Problem 3.a
% define block code
k = 12;
n = 24;
dmin = 8;
% Define Eb/N0 values in dB
%EbN0dB = linspace(-5, 10, 100);
EbN0dB = -12:0.5:12; % Range of Eb/N0 values in dB

% Convert Eb/N0 to linear scale
EbN0 = 10.^((EbN0dB / 10));

% Calculate the theoretical Pe for BPSK
Pb_theoretical = qfunc(sqrt(2*EbN0));

% calculate block error probability using binomial distribution formula

Pe_theoretical = 1 - (1-Pb_theoretical).^n;

% Plot the results
figure()
semilogy(EbN0dB, Pe_theoretical, 'LineWidth', 2);
title('Uncoded Probability of block error Pe vs. Eb/N0, k=12, n=24, dmin=8');
xlabel('Eb/N0 (dB)');
ylabel('P_e');
grid on;
