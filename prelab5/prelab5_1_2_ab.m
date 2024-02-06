%% Prob 2.a

% Define Eb/N0 values in dB
EbN0dB_1 = linspace(-5, 10, 100);

% Convert Eb/N0 to linear scale
EbN0_1 = 10.^(EbN0dB_1 / 10);

% Calculate the theoretical Pe for BPSK
Pe_theoretical = 0.5 * qfunc(sqrt(EbN0_1));
% Plot the results
figure(1);

semilogy(EbN0dB_1, Pe_theoretical, 'LineWidth', 2);
title('Uncoded Probability of Single Data Bit Error vs. Eb/N0 for BPSK');
xlabel('Eb/N0 (dB)');
ylabel('P_e');
grid on;

%% Problem 2.b
% Define Eb/N0 values in dB
EbN0dB_2 = linspace(-5, 10, 100);

% Convert Eb/N0 to linear scale
EbN0_2 = 10.^(EbN0dB_2 / 10);

% Calculate the theoretical Pe for BPSK
Pb_theoretical_2 = 0.5 * qfunc(sqrt(EbN0_2*0.5));
% Plot the results
hold on
semilogy(EbN0dB_2, Pb_theoretical_2, 'LineWidth', 2);
title('Uncoded Probability of Single Data Bit Error vs. Eb/N0 for BPSK');
xlabel('Eb/N0 (dB)');
ylabel('P_e');
grid on;
legend('Rc=1','Rc=1/2')