function thresholded_vector = binary_threshold(input_vector)
    % Perform binary thresholding. If an element is greater than 0 set it
    % equal to 1, else all other elements equal zero.
    thresholded_vector = (input_vector > 0);
end