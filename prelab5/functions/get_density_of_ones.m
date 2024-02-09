function [relative_density] = get_density_of_ones(data)
%GET_DENSITY_OF_ONES Summary of this function goes here
% Count the number of ones and zeros
num_ones  = sum(data == 1);
num_zeros = sum(data == 0);

% Calculate the relative density
relative_density = num_ones / num_zeros;

% Display the results
fprintf('   Number of Ones: %d\n', num_ones);
fprintf('   Number of Zeros: %d\n', num_zeros);
fprintf('   Relative Density (Ones to Zeros): %.2f\n', relative_density);
end

