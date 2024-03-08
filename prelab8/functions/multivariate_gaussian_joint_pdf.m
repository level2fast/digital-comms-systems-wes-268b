function joint_pdf = multivariate_gaussian_joint_pdf(X, mu, sigma)
    n = size(X, 1);  % Number of random variable vectors
    joint_pdf = zeros(n, 1);  % Initialize joint PDF vector
 %   for i = 1:n
        x = X';  % Extract random variable vector
        exponent = -0.5 * (x - mu)' * inv(sigma) * (x - mu);
        coef = 1 / sqrt((2 * pi)^size(x, 1) * det(sigma));
        factor = exp(exponent);
        joint_pdf = coef * factor;
%    end
end