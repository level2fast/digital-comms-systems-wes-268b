%% OFDM and the FFT Matrix
%% 2.a Show that all columns are orthogonal
% Create an IDFT matrix
[twiddle_factors_mat] = calc_dft_twiddle_factors([1 1 1 1]);
disp('IDFT Matrix:');
disp(twiddle_factors_mat);
A = twiddle_factors_mat;
n = length(twiddle_factors_mat);
is_col_orthogonal = false;
% Iterate over pairs of columns
for i = 1:n-1
    for j = i+1:n
        % Calculate the dot product
        dot_product = dot(A(:,i), A(:,j));
        % Check if the dot product is close to zero (considering numerical precision)
        if abs(dot_product) > 1e-10
            is_col_orthogonal = false;
            return;
        end
    end
end
msg = sprintf('2.a: Dot product of Twiddle Factor columns: %d \n', is_col_orthogonal);
disp(msg);
identity = twiddle_factors_mat *conj(transpose(twiddle_factors_mat));
identity = identity/norm(identity);




%% 2.b.i Verify that FF^-1 = I, in other words, show that the matrix F is the inverse of F^-1;
[dft_twiddle_factors_mat] = calc_dft_twiddle_factors([1 1 1 1]);
[idft_twiddle_factors_mat] = calc_idft_twiddle_factors([1 1 1 1]);
disp('2.b.i: DFT Matrix * IDFT Matrix):');
disp(dft_twiddle_factors_mat * idft_twiddle_factors_mat);



