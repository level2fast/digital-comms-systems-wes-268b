%% Problem 2.1.a Simulation of a Scrambler

% load raw data that will be used for scrambler
raw_data = load("RawData.mat");

%% Determine the relative density of logic ones to logic zeros
disp('Determine the relative density of logic ones to logic zeros')
data = vertcat( raw_data.rawdata );

% Count the number of ones and zeros
num_ones  = sum(data == 1);
num_zeros = sum(data == 0);

% Calculate the relative density
relative_density = num_ones / num_zeros;

% Display the results
fprintf('   Number of Ones: %d\n', num_ones);
fprintf('   Number of Zeros: %d\n', num_zeros);
fprintf('   Relative Density (Ones to Zeros): %.2f\n', relative_density);

%% Determine the longest run of Ones and longest run of zeros

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
fprintf('   The longest run of 1''s in the vector is: %d\n\n', max_run);
%% Scramble this data set by using bitwise XOR of this data with a 2^5-1 PRBS. Use pngen()
disp('Scramble this data set by using bitwise XOR of this data with a PRBS');
r = 5;
N = length(data);
prbs = (pngen(r,N));
data_scrambled = xor(data, prbs);
msg  = sprintf("   Data scrambled using pngen(r=%d,N=%d) \n",r,N);
disp(msg)
%% Determine the density of Ones after scrambling and the longest run of Ones and the longest run of Zero's. Explain your result.
disp('Density of Ones after scrambling and the longest run of Ones and the longest run of Zero')
data_density = get_density_of_ones(data_scrambled);
data_lngst_run_ones = get_longest_run_of_ones(data_scrambled);

% The relative density of ones and the longest run of ones increased due to
% randomness of the data 

%% Uncsramble the data by applying an XOR operation with the same PRBS again and confirm that you recover the original data set. 
disp('Unscramble data and recover original signal')
data_unscrambled = xor(data_scrambled, prbs);
if data_unscrambled == data
    disp('  The data was unscrambled.');
else
    disp('  The data is still scrambled.');
end

%% Offset or mis-align the unscrambling block.
disp('Offset or mis-align the unscrambling block.')
% Offset or mis-align the unscrambling block by one bit with respect to the
% scrambling block and compare the new "uncscrambled" block to the original
% data. Explained what happened and why. 
% Example vector
original_vector = data_unscrambled;

% Shifting elements to the right by 1 index
shifted_vector = [original_vector(2:end), 0];
data_unscrambled_offset = shifted_vector; % shifted unscrambled data 1 bit to the left.

if data_unscrambled_offset == data
    disp('  data_unscrambled_offset == data');
else
    disp('  data_unscrambled_offset != data');
end

% Answer: The data values are not equal. Using a scrambler does not gaurantee 
% proper encoding of the data. If data bits become misaligned errors will
% occur.
