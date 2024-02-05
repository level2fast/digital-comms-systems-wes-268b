%% Problem 2.1.a Simulation of a Scrambler

% load raw data that will be used for scrambler
raw_data = load("RawData.mat");

% Deetermine the relative density of logic ones to logic zeros =
data = vertcat( raw_data.rawdata );

% Count the number of ones and zeros
num_ones = sum(data == 1);
num_zeros = sum(data == 0);

% Calculate the relative density
relative_density = num_ones / num_zeros;

% Display the results
disp('Raw Dataset:');
fprintf('Number of Ones: %d\n', num_ones);
fprintf('Number of Zeros: %d\n', num_zeros);
fprintf('Relative Density (Ones to Zeros): %.2f\n', relative_density);

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
fprintf('The longest run of 1''s in the vector is: %d\n', max_run);




% Determine the longest run of Ones and longest run of zeros

% Scramble this data set by using bitwise XOR of this data with
% a 2^5-1 PRBS. Use pngen()

% Determine the density of Ones after scrambling and the longest run of
% Ones and the longest run of Zero's. Explain your result.

% Uncsramble the data by applying an XOR operation with the same PRBS again
% and confirm that you recover the original data set. 

% Offset or mis-align the unscrambling block by one bit with respect to the
% scrambling block and compare the new "uncscrambled" block to the original
% data. Explained what happened and why. 
