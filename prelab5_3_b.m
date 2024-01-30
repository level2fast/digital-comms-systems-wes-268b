%% Problem 3.a
% MATLAB script for plotting Pe vs. Eb/N0 using the binomial distribution formula

% Parameters
Eb_N0_dB = -10:0.5:10; % Range of Eb/N0 values in dB
Eb_N0 = 10.^(Eb_N0_dB / 10); % Convert dB to linear scale

% Transmission parameters
n = 24; % Number of bits in a block
k = 12; % Number of information bits
p = 0.1; % Probability of bit error

% Calculate Pe using binomial distribution formula
Pe = zeros(size(Eb_N0));
for i = 1:length(Eb_N0)
    % Calculate noise variance (sigma^2) from Eb/N0
    sigma_square = Eb_N0(i) / (2 * k);
    
    % Calculate Pe using the binomial distribution formula
    Pe(i) = 1 - binocdf(floor((n-k)/2), n, p);
end

% Plot Pe vs. Eb/N0
figure;
semilogy(Eb_N0_dB, Pe, 'b-o', 'LineWidth', 2);
grid on;
title('Probability of Block Error vs. Eb/N0');
xlabel('Eb/N0 (dB)');
ylabel('Probability of Block Error (Pe)');
legend('Binomial Distribution');
