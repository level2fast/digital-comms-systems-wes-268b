% MATLAB script for plotting Pe vs. Eb/N0 using the binomial distribution formula
% Parameters
Eb_N0_dB = -10:0.5:10; % Range of Eb/N0 values in dB
Eb_N0 = 10.^(Eb_N0_dB / 10); % Convert dB to linear scale

% Transmission parameters
n = 24; % Number of bits in a block
k = 12; % Number of information bits
R_c = k/n;

Pe = zeros(size(Eb_N0));

% Calculate the theoretical probability of a single bit erorr for BPSK
Pb_theoretical = qfunc(sqrt(2*Eb_N0*R_c));

% The probability of block error can be found using a cumulative binomial 
% distribution formula. This loop calculates the probability of of block
% error at a specific SNR 
for idx = 1:length(Eb_N0)
    % Calculate the pobility of a data block having 12 bits correct given the
    % probability of a single data bit eror for BPSK. This will tell us the
    % probability of a block error given that we have n number of bits in a 
    % block(aka trials) and k information bits(i.e. expected possible successes).
    Pe(idx) = 1 - binocdf(floor((n-k)/2), n, Pb_theoretical(idx));
end

% Plot Pe vs. Eb/N0
figure(1);
semilogy(Eb_N0_dB, Pe, 'b-o', 'LineWidth', 2);
grid on;
title('Probability of Block Error vs. Eb/N0');
xlabel('Eb/N0 (dB)');
ylabel('Probability of Block Error (Pe)');
legend('Binomial Distribution');
