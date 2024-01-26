function [vectors_16QAM_gray_code] = create_16_qam_gray_code(print)
%CREATE_16_QAM_GRAY_CODE Summary of this function goes here
Define the 16-QAM Gray code mapping
gray_code_mapping = [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8];

Convert each Gray code index to 4-bit binary representation
binary_gray_code = dec2bin(gray_code_mapping, 4) - '0';

Create 16 MATLAB vectors representing 16-QAM Gray code
vectors_16QAM_gray_code = zeros(16, 4);

for i = 1:16
    vectors_16QAM_gray_code(i, :) = binary_gray_code(i, :);
end

if(print ==0 )
    disp('16-QAM Gray code vectors:');
    disp(vectors_16QAM_gray_code);
end
end

