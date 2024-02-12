function isOrthogonal = is_orthogonal(A)
    % Check if matrix A is symmetric
    is_conj_symmetric = isequal(A,conj(transpose(A)));
    if ~is_conj_symmetric
        error('Matrix is not symmetric.');
    end

    % Get the size of the matrix
    [n, m] = size(A);

    % Ensure the matrix is square
    if n ~= m
        error('Matrix is not square.');
    end
    
    % Initialize flag for orthogonality
    isOrthogonal = true;

    % Iterate over pairs of columns
    for i = 1:n-1
        for j = i+1:n
            % Calculate the dot product
            dot_product = dot(A(:,i), A(:,j));

            % Check if the dot product is close to zero (considering numerical precision)
            if abs(dot_product) > 1e-10
                isOrthogonal = false;
                return;
            end
        end
    end
end
