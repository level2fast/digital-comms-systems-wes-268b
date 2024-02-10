% This script produces a simulation of a QPSK modulator/demodulator with coding. 

num_symbols = 16;  % Number of QPSK symbols
rng(42);  % Set the random seed for reproducibility
k = 1;    % information bits
n = 3;    % total bits
Rc = k/n; % code rate

% This assumes that the transmitted pulse shape is rectangular and 
% modulation mapping for each binary data stream
Eb_N0_dB = -1:1:12;                % Range of Eb/N0 values in dB
Eb_N0    = 10.^(Eb_N0_dB / 10);    % Convert dB to linear scale
BER      = zeros(1,length(Eb_N0)); % pre-allocate BER vector
N        = 1;                      % samples per symbol

i_decoded_data = zeros(1,num_symbols*N/2);    % data vector 
q_decoded_data = zeros(1,num_symbols*N/2);    % data vector 
sigma = 1/sqrt(2);                            % noise variance

% for each signal to noise ratio
for idx=1:length(Eb_N0)
    num_errors = 0;
    num_bits   = 0;
    symbol_error = 0;
    i_error = 0;
    q_error = 0;
    % calc amplitude of qpsk symbol
    A = sqrt(Rc * Eb_N0(idx)); 
    while( num_errors < 200 )
        %% Modulator
        % Create  Data source which consists of a random sequence of 0's and 1's
        % Parameters
        data_bits = randi([0, 1], 1, 2 * num_symbols);  % Random binary data
        
        % QPSK Modulation
        modulated_symbols = qpsk_modulate(data_bits); 
       
        %  Separate the data source into two data streams for the I and Q channels
        % Separate data into I and Q channels
        I_channel = real(modulated_symbols(1:2:end));
        Q_channel = real(modulated_symbols(2:2:end));
        
        i_rect_pulse_uncoded = I_channel;
        i_rect_pulse_uncoded(i_rect_pulse_uncoded > 0) = 1;
        i_rect_pulse_uncoded(i_rect_pulse_uncoded < 0) = 0;

        q_rect_pulse_uncoded = Q_channel;
        q_rect_pulse_uncoded(q_rect_pulse_uncoded > 0) = 1;
        q_rect_pulse_uncoded(q_rect_pulse_uncoded < 0) = 0;

        % Channel code each stream using the repeat (3,1) code
        i_modulated_symbols_encoded = repetition_code_encode(i_rect_pulse_uncoded);
        q_modulated_symbols_encoded = repetition_code_encode(q_rect_pulse_uncoded);
        
        % Scramble the coded sequence using a scrambler
        [i_scrambled_data, i_prbs] = scramble_data(i_modulated_symbols_encoded,5);
        [q_scrambled_data, q_prbs] = scramble_data(q_modulated_symbols_encoded,5);

        % Threshold data and apply pulse shaping
        i_rect_pulse = i_scrambled_data;
        i_rect_pulse = (i_rect_pulse * 2) - 1;
        i_rect_pulse = i_rect_pulse * A;

        q_rect_pulse = q_scrambled_data;
        q_rect_pulse = (q_rect_pulse * 2) -1;
        q_rect_pulse =  q_rect_pulse * A;

        %% Channel 
        % Generate independent zero-mean Gaussian random variable with variance of 1/2
        noise = sigma * randn(1, length(i_rect_pulse));
        i_rect_pulse_noisy = noise + i_rect_pulse;
        q_rect_pulse_noisy = noise + q_rect_pulse;

        %% Receiver
        % Decode each of the unscrambled binary streams separately using a
        % majority logic rule.
        i_rect_pulse_decoded = i_rect_pulse_noisy;
        i_rect_pulse_decoded(i_rect_pulse_decoded > 0) = 1;
        i_rect_pulse_decoded(i_rect_pulse_decoded < 0) = 0;

        q_rect_pulse_decoded = q_rect_pulse_noisy;
        q_rect_pulse_decoded(q_rect_pulse_decoded > 0) = 1;
        q_rect_pulse_decoded(q_rect_pulse_decoded < 0) = 0;

        % unscramble the i/q data
        i_data_unscrambled = xor(i_rect_pulse_decoded, i_prbs);
        q_data_unscrambled = xor(q_rect_pulse_decoded, q_prbs);

        % Now use majority logic rule to check every 3 bits. 
        % loop through i data 
        len = length(i_data_unscrambled)/3; % determine how many steps we should take for every 3 bits
        end_idx = 0;
        start_idx = 1;
        i_decoded_data(:) = 0;
        q_decoded_data(:) = 0;
        for idxDecode=1:len
            end_idx = 3* idxDecode; % calculate first end index

            i_three_bits = i_data_unscrambled(start_idx:end_idx);
            i_decoded_data(idxDecode) = majority_logic_rule(i_three_bits);
            
            q_three_bits = q_data_unscrambled(start_idx:end_idx);
            q_decoded_data(idxDecode) = majority_logic_rule(q_three_bits);
            
            start_idx = start_idx + 3;
        end
        
        % Perform element-wise comparison of bits and sum the number of bits that are not equal
        i_error = compare_and_sum(i_decoded_data,i_rect_pulse_uncoded);
        q_error = compare_and_sum(q_decoded_data,q_rect_pulse_uncoded);

        % compute the total number of errors
        num_errors = num_errors + i_error + q_error;

        % compute the probability of symbol error
        % look at every 3 bits of the decoded message and determine the
        % symbol error rate
        for idxSymblErr=1:length(i_decoded_data)
            % extract 1 symbol from decoded vector i.e. i and q component
            i_comp_dec = i_decoded_data(idxSymblErr);
            q_comp_dec = q_decoded_data(idxSymblErr);
            symbol_decoded = [i_comp_dec q_comp_dec];

            % extract 1 symbol from uncoded vector i.e. i and q component
            i_comp_unc = i_rect_pulse_uncoded(idxSymblErr);
            q_comp_unc = q_rect_pulse_uncoded(idxSymblErr);
            symbol_uncoded = [i_comp_unc q_comp_unc];

            % compare the symbols
            if(isequal(symbol_decoded,symbol_uncoded))
                symbol_error = symbol_error + 1;
            end
        end
        % accumulate errors / bits
        num_bits = num_bits + num_symbols;
    end 
   BER(idx) = num_errors/num_bits;
end

%% Q2.2.1 - Plot the uncoded probability of a bit error down to an error rate of about 10e−4
% MATLAB script for plotting Pe vs. Eb/N0 for QPSK
% QPSK modulation parameters
M = 4;       % QPSK modulation order
k = log2(M); % Number of bits per symbol

% Calculate Pe for QPSK
Pe_uncoded = qfunc(sqrt(2*Eb_N0 / k)); % Uncoded BER for QPSK

% Plot Pe vs. Eb/N0
figure(1);
semilogy(Eb_N0_dB, Pe_uncoded, 'b-o', 'LineWidth', 2);
grid on;
title('Uncoded Probability of Single Data Bit Error (Pe) vs. Eb/N0 for QPSK');
xlabel('Eb/N0 (dB)');
ylabel('Probability of Single Data Bit Error (Pe)');
legend('Uncoded QPSK');
grid on;

%% Q2.2.1 - Plot the uncoded probability of a symbol error for QPSK down to an error rate of about 10e−4
% MATLAB script for plotting Pe vs. Eb/N0 for uncoded QPSK
% QPSK modulation parameters
M = 4;       % QPSK modulation order
k = log2(M); % Number of bits per symbol

% Calculate Pe for QPSK (uncoded symbol error rate)
energy_per_symbol = 2*qfunc(sqrt(2*Eb_N0 / k)); % energy per symbol is 2 * energy per bit for qpsk
Pe_symbol_uncoded = 1 - (1 - energy_per_symbol);

% Plot Pe vs. Eb/N0
figure(2);
semilogy(Eb_N0_dB, Pe_symbol_uncoded, 'b-o', 'LineWidth', 2);
grid on;
title('Uncoded Probability of Symbol Error (Pe) vs. Eb/N0 for QPSK');
xlabel('Eb/N0 (dB)');
ylabel('Probability of Symbol Error (Pe)');
legend('Uncoded QPSK');

%% Q2.2.2 Plot the coded probability of a bit error for one quadrature components using the (3, 1) repeat code.
figure(3)
semilogy(Eb_N0_dB,BER,'b-o','LineWidth',2);
grid on;
xlabel('EbNo (dB)');
ylabel('BER');
title('BER vs. EbNo (dB)');
legend('Simulated BER');


%% Q2.2.3 Determine the coding gain at pe = 10−4 and compare this simulated value with theory for BPSK for a single quadrature component.
figure(4);
semilogy( Eb_N0_dB,BER,'b-o',Eb_N0_dB, Pe_uncoded,'r-o',"LineWidth",2);
grid on;
xlabel('EbNo (dB)');
ylabel('BER');
title('BER vs. EbNo (dB)');
legend('Simulated BER','Theory BER');

Display("Q2.2.3 Answer: Coding gain is approximately 2.5dB")
