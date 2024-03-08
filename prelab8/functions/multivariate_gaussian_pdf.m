function pdf_value = multivariate_gaussian_pdf(x, mu, sigma)
    n = length(x);
    exponent = -0.5 * (x - mu).' * inv(sigma) * (x - mu);
    coef = 1 / (sqrt((2 * pi)^n * det(sigma)));
    pdf_value = coef * exp(exponent);
end


