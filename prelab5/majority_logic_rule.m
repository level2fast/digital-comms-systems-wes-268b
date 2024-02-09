function result = majority_logic_rule(data)
    % Input:
    %   data: Input block of binary data (row vector of 0s and 1s)
    % Output:
    %   result: Result of the majority logic rule (0 or 1)

    % Check if the input data is a row vector
    arguments
        data (1,3) {mustBeNonempty}
    end
    if ~isrow(data)
        error('Input data must be a row vector.');
    end
    bit_0 = data(1);
    bit_1 = data(2);
    bit_2 = data(3);

    first_and_gate  = ~bitand(bit_0,bit_1);
    second_and_gate = ~bitand(bit_1,bit_2);
    third_and_gate  = ~bitand(bit_0,bit_2);
    fourth_and_gate = bitand(first_and_gate,third_and_gate);
    final_and_gate  = ~bitand(second_and_gate,fourth_and_gate);
    % Perform majority logic rule
    result = final_and_gate;

    % disp('Input Data:');
    % disp(data);
    % disp(['Result of Majority Logic Rule: ' num2str(result)]);
end
