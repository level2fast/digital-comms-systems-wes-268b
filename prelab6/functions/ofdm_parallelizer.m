function [complex_block] = ofdm_parallelizer(ofdm)
%OFDM_PARALLELIZER Summary of this function goes here
%   Detailed explanation goes here
arguments
    ofdm.symbol_values  (1,:) {mustBeNonempty} = 0
    ofdm.subcarrier_num (1,1) {mustBeNonempty} = 32
end
symbols = ofdm.symbol_values;
subcarriers = ofdm.subcarrier_num;
% Reshape the vector into a matrix with 32 rows
complex_block = reshape(symbols, subcarriers, []);

% Display the reshaped matrix
disp('Original vector:');
disp(symbols);
disp('Reshaped matrix with 32 rows:');
disp(complex_block);
end

