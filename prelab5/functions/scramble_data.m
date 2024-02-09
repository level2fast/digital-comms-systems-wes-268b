function [data_scrambled, prbs] = scramble_data(data,pn_order)
%SCRAMBLE_DATA Summary of this function goes here
%   Detailed explanation goes here
r = pn_order;
N = length(data);
prbs = (pngen(r,N));
data_scrambled = xor(data, prbs);
end

