
%%1.e Union bound problem

% calculate the minimum hamming weight of all code words in a 16-QAM
% binary vector sequence
codewords = [...
     0     0     0     0;
     0     0     0     1;
     0     0     1     1;
     0     0     1     0;
     0     1     1     0;
     0     1     1     1;
     0     1     0     1;
     0     1     0     0;
     1     1     0     0;
     1     1     0     1;
     1     1     1     1;
     1     1     1     0;
     1     0     1     0;
     1     0     1     1;
     1     0     0     1;
     1     0     0     0;];

 dmin = 0;
 temp = 0;
 for i = 1:length(codewords)
     for j = 1:length(codewords)
         % Try to determine the lowest hamming distance between all code
         % words in our 16-QAM constellation
        temp = calculate_min_hamming_distance(codewords(i),codewords(j));
        if(temp > 0)
            dmin = temp;
        end
     end
     temp = 0;
 end
msg = sprintf("16-QAM dmin:%i",dmin);
display(msg)

