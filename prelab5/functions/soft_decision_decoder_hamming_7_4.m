function [decoded_message] = soft_decision_decoder_hamming_7_4(data,H,k_arg,n_arg)
%SOFT_DECISION_DECODER Summary of this function goes here
% Takes the received codeword and correlate against every possible codeword
% in the codeword space and pick the one that has the best correlation. 
% Computes every possible code word in the codebook and construct it as matrix 

recv_codeword = data; % received coded sequence
k = k_arg; % number of information bits
n = n_arg;
%% compute the correlation vector
G = get_generator_matrix(H); % generator matrix
valid_codes_mat = zeros(2^k-1,n); % matrix of valid code words sequence dimensions [2^k x n]
for validIdx = 0:2^k-1
	m_vec = base2dec(dec2bin(validIdx,k).',2).';
	valid_codes_mat(validIdx+1,:) = mod(m_vec * G,2);
end

% Soft decision decoder
corr_val = recv_codeword * (2 * valid_codes_mat.'- 1);
[val, idx]  = max(corr_val,[],2); % from 2^k correlation values, the index 
                                  % where the correlation value is maximized 
                                  % corresponds to the maximum likelihood 
                                  % transmit code word.

decoded_message = valid_codes_mat(idx,1:k);
end

