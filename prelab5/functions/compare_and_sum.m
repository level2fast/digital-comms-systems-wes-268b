function sum_diff_elements = compare_and_sum(vector1, vector2)
    % Input:
    %   vector1: First input vector
    %   vector2: Second input vector
    % Output:
    %   sum_diff_elements: Sum of the number of elements that are not equal

    % Check if the input vectors are of the same size
    if numel(vector1) ~= numel(vector2)
        error('Input vectors must be of the same size.');
    end

    % Perform element-wise comparison and sum the number of elements that are not equal
    sum_diff_elements = sum(vector1 ~= vector2);

    % disp('Vector 1:');
    % disp(vector1);
    % disp('Vector 2:');
    % disp(vector2);
    % disp(['Sum of Different Elements: ' num2str(sum_diff_elements)]);
end
