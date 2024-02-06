function isValid = check_code_word(codeword, H)
    % Check if the codeword is a valid codeword using the parity check matrix H
    syndrome = mod(codeword * H', 2);  % Calculate the syndrome

    % If the syndrome is all zeros, the codeword is valid
    isValid = all(syndrome == 0);
end
