% autocovariance matrix C of a multivariate gaussian probability density function 
C =  reshape([1 1 1 4],2,2);

%% a.) find a matrix R that diagonalizes C
% R is constructed from the normalized eigenvectors of C
[R_matrix,Cd] = eig(C);

% calculate c prime which is a dian
c_prime_matrix = R_matrix*C*((R_matrix.'));
%% b.) Determine an expression for the new joint probability distribution fx′y′(x′, y′) in terms of σx′ and σy′ .
mu = [0 0];
Sigma = C;

x1 = -3:0.2:3;
x2 = -3:0.2:3;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mu,Sigma);
y = reshape(y,length(x2),length(x1));

figure(1)
surf(x1,x2,y)
axis([-3 3 -3 3 0 0.4])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')

xprime_yprime = R_matrix * X.';
mu = 0;
y_decorr = mvnpdf(xprime_yprime',mu,Sigma);
y_decorr = reshape(y_decorr,length(x2),length(x1));
figure(2)
surf(x1,x2,y_decorr)
axis([-3 3 -3 3 0 0.4])
xlabel('x1 prime')
ylabel('x2 prime')
zlabel('Probability Density')
%% c.) Plot and compare the distributions of fxy(x, y) and fx′y′(x′, y′).