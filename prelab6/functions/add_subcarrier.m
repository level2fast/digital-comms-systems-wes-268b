function [output_signal] = add_subcarrier(signal,f_delta,time_vec)
%ADD_SUBCARRIER Summary of this function goes here
%   Detailed explanation goes here
output_signal = signal* exp(1j * 2 * pi * f_delta * time_vec);
end

