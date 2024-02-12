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
n = 4; % You can change this to your desired size

% Create a diagonal matrix with either 1 or -1 on the diagonal
D = diag(randi([0, 1], n, 1) * 2 - 1);

disp('Orthogonal Symmetric Matrix:');
disp(D);

result = is_orthogonal(D);
disp(['The columns of the matrix are orthogonal: ', mat2str(result)]);

