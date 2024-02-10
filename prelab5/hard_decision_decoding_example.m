%% 2.3 Comparison of Hard and Soft Decoding for the Hamming Code

% Load data
% data = load("msg_coded.mat");
% rx_data_noisy = data.rx_blocks;
% H     = data.H;
% x_dim = data.x;
% y_dim = data.y;
% Rc    = data.Rc;
% n     = data.n;
% k     = data.k;

%% Use hard decision detection to decode the message
% Define a (7,4) Hamming Code
% G = [1 1 0 1 0 0 0; 0 1 1 0 1 0 0; 1 1 1 0 0 1 0; 1 0 1 0 0 0 1];
H = [1,0,1,1,1,0,0;
     1,1,0,1,0,1,0;
     0,1,1,1,0,0,1];

[row, col] = size(H);

% Calculate the systematic generator matrix G
k = col - row;  % Number of information bits
I_k = eye(k);   % Form the identity matrix for k bits

% Create a systematic generator matrix G
K_t = H(:, 1:k)'; % grab information bits from parity check matrix
G = [I_k K_t];    % combine identity matrix with information bits transposed 

% Original message
message = [1 0 1 0];

% Encoding: Multiply the message by the generator matrix
codeword = mod(message * G, 2);

% Add some simulated errors to the codeword
received_codeword = codeword;
received_codeword(3) = mod(received_codeword(3) + 1, 2);
decoded_message = hamming_7_4_decoder(received_codeword,H);

% Display results
disp('Original Message:');
disp(message);
disp('Encoded Codeword:');
disp(codeword);
disp('Received Codeword:');
disp(received_codeword);
disp('Decoded Codeword:');
disp(decoded_message);
fprintf("Most likely transmitted message\n");
fprintf('decoded message vector: [')
fprintf('%d, ', decoded_message(1:3))
fprintf('%d]\n\n', decoded_message(4))
