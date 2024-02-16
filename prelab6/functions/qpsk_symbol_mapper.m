function [comp_baseband_qpsk_sym] = qpsk_symbol_mapper(data)
%SYMBOL_MAPPER Takes 2 bits and creates one complex baseband QPSK symbol.
% Map binary data to QPSK symbols
arguments
    data.two_bits (1,2) {mustBeInRange(data.two_bits,0,1)} 
    data.bit_vector (1,:) {mustBeInRange(data.bit_vector,0,1)}
    data.input {mustBeMember(data.input,{'bit','bitstream'})} = 'bitstream'
end
if(data.input == "bit")
    qpsk_symbols = 1/sqrt(2) * (2 * bit1 - 1 + 1j * (2 * bit2 - 1));
    comp_baseband_qpsk_sym = qpsk_symbols;
else
    % Map binary data to QPSK symbols
    comp_baseband_qpsk_sym = 1/sqrt(2) * (2 * data.bit_vector(1:2:end) - 1 + 1j * (2 * data.bit_vector(2:2:end) - 1));

end

