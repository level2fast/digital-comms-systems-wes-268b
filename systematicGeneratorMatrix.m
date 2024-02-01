function G = systematicGeneratorMatrix(H)
    % Check if H is a valid parity matrix
    [m, n] = size(H);
    if m >= n
        error('Invalid parity matrix. Number of rows should be less than the number of columns.');
    end
    
    % Calculate the systematic generator matrix G
    k = n - m;  % Number of information bits
    I_k = eye(k);  % Identity matrix of size k
    
    % Create a systematic generator matrix G
    G = [I_k H(:, 1:k)'];
end
