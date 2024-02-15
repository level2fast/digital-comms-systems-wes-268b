function [comp_baseband_qpsk_sym] = symbol_mapper(bit1,bit2)
%SYMBOL_MAPPER Takes 2 bits and creates one complex baseband QPSK symbol.
% Map binary data to QPSK symbols
qpsk_symbols = 1/sqrt(2) * (2 * bit1 - 1 + 1j * (2 * bit2 - 1));
comp_baseband_qpsk_sym = qpsk_symbols;
end

