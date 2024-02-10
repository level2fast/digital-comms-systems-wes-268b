function decoded_message = hamming_7_4_decoder(encoded_message,parity_check_mat)
    % Function to decode a message block using Hamming (7,4) code. This
    % function uses hard decision decoding to decode the incoming message.
    % In hard decision decoding for a Hamming(7,4) code, you make decisions 
    % based on the received bits without considering the reliability or 
    % likelihood of each bit. The basic idea is to use the received bits 
    % directly and apply error correction if needed.

    % Check if the input length is a multiple of 7
    if mod(length(encoded_message), 7) ~= 0
        error('Input length must be a multiple of 7');
    end
    
    encoded_message_noisy = encoded_message;
    H = parity_check_mat;

    % Check if the codeword is a valid codeword using the parity check matrix H
    syndrome = mod(encoded_message_noisy * H', 2);  % Calculate the syndrome
    % fprintf('syndrome vector: [')
    % fprintf('%d, ', syndrome(1:end-1))
    % fprintf('%d]\n', syndrome(end))

    % If the syndrome is all zeros, the codeword is valid
    isValid = all(syndrome == 0);
    
    if isValid
        %disp('The codeword is valid.');
        decoded_message = encoded_message_noisy(1:4);
        return
    else
        %disp('The codeword is invalid.');
    end

    % Perform error correction
    %errorPosition = bin2dec(num2str(flip(syndrome)));
    error_position = (length(encoded_message_noisy) - bin2dec(num2str(flip(syndrome)))) +1;

    % Correct the error, TODO sdd: I need to index starting from the right
    % most position in this vector, matlab defaults to left most
    encoded_message_noisy(error_position) = mod(encoded_message_noisy(error_position) + 1, 2);

    % Extract 1st 4 bits which correspond to the message that was sent 
    decoded_message = encoded_message_noisy(1:4);
    
    % fprintf("Most likely transmitted codeword\n");
    % fprintf('decoded codeword vector: [')
    % fprintf('%d, ', decoded_message(1:end-1))
    % fprintf('%d]\n\n', decoded_message(end))
end