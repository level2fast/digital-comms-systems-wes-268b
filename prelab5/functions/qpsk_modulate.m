% QPSK Modulation Function
function modulated_symbols = qpsk_modulate(data_bits)
    % QPSK modulation mapping
    constellation = exp(1j * pi/4 * [0 1 2 3]);

    % Reshape data bits into symbol pairs
    symbol_pairs = reshape(data_bits, 2, length(data_bits)/2);

    % Modulate symbol pairs using QPSK mapping
    symbol_pairs_transposed = symbol_pairs.';
    symbol_pairs_mapped = bi2de(symbol_pairs_transposed) + 1;
    modulated_symbols = constellation(symbol_pairs_mapped);
end