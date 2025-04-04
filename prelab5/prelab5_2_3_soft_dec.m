%% 2.3 Soft Decoding for the Hamming Code
% Load data
data = load("msg_coded.mat");
received_codeword = data.rx_blocks;
H     = data.H;
x_dim = data.x;
y_dim = data.y;
Rc    = data.Rc;
n     = data.n;
k     = data.k;

%% Use hard decision detection to decode the message
% First apply thresholding to received codewords
received_codeword = binary_threshold(received_codeword);
decoded_message = zeros(size(received_codeword,1),k);
iterations = length(received_codeword);
for idx=1:iterations
    codeword_row = received_codeword(idx,:);
    decoded_message(idx,:) = soft_decision_decoder_hamming_7_4(codeword_row,H,k,n);
end

% reformat matrix to have 8 columns 1 for each bit in the image
new_mat = zeros(size(received_codeword,1)/2,k*2);
cur_row  = 1;
next_row = 2;
iter_count = (size(decoded_message,1)/2);
temp = decoded_message;
for idxNewMat = 1:iter_count
    row1 = temp(cur_row,:);
    row2 = temp(next_row,:);
    
    % Add combined rows to new matrix
    new_mat(idxNewMat,:) = [row1, row2];

    cur_row = cur_row + 2;
    next_row = next_row + 2;
end

dec_col_vec = uint8(zeros(size(decoded_message,1)/2,1));
for idxDecMat = 1:iter_count
    dec_col_vec(idxDecMat) = uint8(bi2de(new_mat(idxDecMat,:)));
end

image = (reshape(dec_col_vec,x_dim,y_dim));
imshow(image);

% 6.Now, you will compare the differences between using a soft-decision decoder vs. a hard
% decision decoder
% Answer: 
% The term "hard decision" refers to the discrete and deterministic
% nature of the decoding process, where each received symbol is decisively 
% classified as one of the possible transmitted symbols. This is in contrast
% to "soft decision" decoding, where the decoder considers the reliability 
% or likelihood of each received symbol, often represented as probabilities.
% Softe decision decoding provides more precision at the cost of
% complexity.