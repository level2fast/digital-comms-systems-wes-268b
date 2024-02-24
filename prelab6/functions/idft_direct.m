function [idft] = idft_direct(signal)
% Function to perform Inverse Discrete Fourier Transform
% X - Input frequency domain sequence (vector)
% x - Output time domain sequence (vector)

N = length(signal); % Number of points in the sequence
x = zeros(1, N); % Initialize the output sequence

for n = 0:N-1
    for k = 0:N-1
        x(n+1) = x(n+1) + signal(k+1) * exp(1j * 2 * pi * k * n / N);
    end
end

idft = x / sqrt(N); % Normalize by the number of points
end

