function encoded_data = repetition_code_encode(data)
    % This function encodes a block of data using a (3,1) repitition code.
    %
    % A (3,1) repetition code  is a simple error-correcting code in coding 
    % theory. In this code, each bit of the original message is repeated 
    % three times. The purpose of repetition is to introduce redundancy in 
    % the transmission, allowing for the detection and correction of errors.
    %
    % Input:
    %   data: Input block of data (row vector of 0s and 1s)
    % Output:
    %   encoded_data: Encoded data using (3,1) repetition code

    % Check if the input data is a row vector
    if ~isrow(data)
        error('Input data must be a row vector.');
    end

    % Replicate each bit three times
    encoded_data = repmat(data, 3, 1);

    % Reshape to make it a row vector
    encoded_data = reshape(encoded_data, 1, []);

    % disp('Original Data:');
    % disp(data);
    % disp('Encoded Data:');
    % disp(encoded_data);
end
