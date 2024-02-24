function [ofdm_serial_vector] = ofdm_serializer(ofdm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments
    ofdm.symbol_subcarrier_mat (:,:) {mustBeNonempty} = 0
end
ofmd_mat = ofdm.symbol_subcarrier_mat;

% Reshape the matrix into a vector representing the serial bits
ofdm_serial_vector = reshape(ofmd_mat, 1, []);
end

