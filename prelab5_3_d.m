% calculate the asymptotic coding gain 
% Transmission parameters
n = 24; % Number of bits in a block
k = 12; % Number of information bits
R_c = k/n;
dmin=8;
t = floor(dmin-1/2);
msg1 = sprintf("Maximum number of errors code can correct (t): %i",t);
display(msg1)
G = R_c * (t + 1);

msg = sprintf("Asymptotic coding Gain: %i",G);
display(msg)