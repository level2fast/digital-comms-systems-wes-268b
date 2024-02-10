%% Problem 1.4.a
% MATLAB script for constructing a systematic Generator matrix G. 

% What is a Generator matrix?
%   A generator matrix is a matrix whose rows form a basis for a linear code.
%   It is comprised of an identity matrix combined with an arbitrary matrix. 
%   Multiplying a message in row matrix form by a generator matrix produces a
%   codeword.
%   G = [I_k | A]  where (I) is the identity matrix of k bits and (A) is an
%   arbitrary matrix

% What is a parity check matrix?
%   A parity-check matrix of a linear block code C is a matrix which 
%   describes the linear relations that the components of a codeword must 
%   satisfy. It can be used to decide whether a particular vector is a
%   codeword.
%   H = [transpose(A) | I_n-k], where (A) is the arbitrary matrix and the (I)
%   is the identity matrix composed of (n-k) bits

% Lets say we are given a parity check matrix (H). 
H = [1,0,1,1,1,0,0;
     1,1,0,1,0,1,0;
     0,1,1,1,0,0,1];
% How can we obtain the generator matrix (G)?

% The generator matrix can be obtained from the parity check matrix (H) by
% determining the number of information bits in (H) and then producing an
% identity matrix of "k" information bits (I_k). That identity matrix can
% then be combined with the transpose (rows,cols) in (H) that hold the
% information bits. Reecall (H) is made up of an arbitrary matrix which
% holds codewords (A) and an Identity matrix (I_n-k).

% Check if H is a valid parity matrix, number of rows can not be greater
% than or equal to number of columns
[row, col] = size(H);
if row >= col
    error('Invalid parity matrix. Number of rows should be less than the number of columns.');
end

% Calculate the systematic generator matrix G
k = col - row;  % Number of information bits
I_k = eye(k);   % Form the identity matrix for k bits

% Create a systematic generator matrix G
K_t = H(:, 1:k)'; % grab information bits from parity check matrix
G = [I_k K_t];    % combine identity matrix with information bits transposed 
disp('Systematic Generator Matrix G:');
disp(G);