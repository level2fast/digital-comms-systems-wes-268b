function [IDft] = dft_direct(signal)
%UNTITLED Summary of this function goes here
%   Peforms direct computation of IDFT

N = length(signal);
Q = 2*pi/N;
x_dft = zeros(N);
W = zeros(N,N);
for k = 1:N
    S = 0;
    for n = 1:N
        W(k,n) = exp(-1j * Q * (k-1) * (n-1));
        S = S + W(k,n) * signal(n);
    end
    x_dft(k) = S;
end
IDft = x_dft;
end

