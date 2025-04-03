function min_distance = calculate_min_hamming_distance(codeword1, codeword2)
% Example:
%   codeword1 = [0 1 1 0];
%   codeword2 = [1 1 1 0];

    % Ensure codewords have the same length
    if length(codeword1) ~= length(codeword2)
        error('Codewords must have the same length.');
    end
    
    % Calculate Hamming distance
    distance = sum(codeword1 ~= codeword2);
    
    % Minimum Hamming distance
    min_distance = distance;
end
