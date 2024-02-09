% QPSK Modulation Function
function modulated_symbols = qpsk_modulate(data_bits)
    % QPSK modulation mapping
    constellation = exp(1j * pi/4 * [0 1 2 3]);

    % Reshape data bits into symbol pairs
    symbol_pairs = reshape(data_bits, 2, length(data_bits)/2);

    % Modulate symbol pairs using QPSK mapping
    modulated_symbols = constellation(bi2de(symbol_pairs.') + 1);
end