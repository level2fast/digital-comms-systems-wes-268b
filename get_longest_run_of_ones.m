function [max_run] = get_longest_run_of_ones(data)
%GET_DENSITY_OF_ONES Summary of this function goes here
% Initialize variables
current_run = 0;
max_run = 0;

% Iterate through the vector
for i = 1:length(data)
    if data(i) == 1
        % Increment the current run length if the element is 1
        current_run = current_run + 1;
    else
        % Update the maximum run length if the current run is longer
        max_run = max(max_run, current_run);
        % Reset the current run length for the next run
        current_run = 0;
    end
end

% Check one last time after the loop ends
max_run = max(max_run, current_run);

% Display the result
fprintf('   The longest run of 1''s in the vector is: %d\n', max_run);
end

