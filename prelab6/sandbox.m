% A = [1 -1; 
%     -1  2;]; % Example symmetric matrix
% result = is_orthogonal(A);
% disp(['The columns of the matrix are orthogonal: ', mat2str(result)]);

% % Generate a random matrix
% A = rand(4); % You can choose the size of the matrix
% 
% % Perform QR decomposition
% [Q, R] = qr(A);
% 
% % Q is an orthogonal matrix
% disp('Orthogonal Matrix Q:');
% disp(Q);
% 
% % Verify orthogonality
% disp('Q'' * Q (Should be close to the identity matrix):');
% disp(Q' * Q);


% Define the size of the matrix
% n = 4; % You can change this to your desired size
% 
% % Create a diagonal matrix with either 1 or -1 on the diagonal
% D = diag(randi([0, 1], n, 1) * 2 - 1);
% 
% disp('Orthogonal Symmetric Matrix:');
% disp(D);
% 
% result = is_orthogonal(D);
% disp(['The columns of the matrix are orthogonal: ', mat2str(result)]);



% Define parameters for the cosine waves
A1 = 1; % Amplitude of the first cosine wave
f1 = 2; % Frequency of the first cosine wave (Hz)
A2 = 1; % Amplitude of the second cosine wave
f2 = 3; % Frequency of the second cosine wave (Hz)
phi1 = 0; % Phase of the first cosine wave (radians)
phi2 = 0; % Phase of the second cosine wave (radians)

% Define the time interval for integration
t_start = 0; % Start time
t_end = 1; % End time

% Define the function to integrate
func = @(t,Fc,Fd) cos(2*pi*Fc*t) + 2*cos(2*pi*(Fc + Fd)*t);

% Perform the integration
integral_value = integral(@(t) func(t,1,1), t_start, t_end);

% Display the result
disp(['Integral value: ', num2str(integral_value)]);

