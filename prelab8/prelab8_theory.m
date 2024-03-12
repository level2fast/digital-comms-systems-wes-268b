% autocovariance matrix C of a multivariate gaussian probability density function 
C =  reshape([1 1 1 4],2,2);

%% a.) find a matrix R that diagonalizes C
% R is constructed from the normalized eigenvectors of C
[R_matrix,Cd] = eig(C);
disp(R_matrix)

% calculate c prime which is a dian
c_prime_matrix = R_matrix*C*((R_matrix.'));
%% b.) Determine an expression for the new joint probability distribution fx′y′(x′, y′) in terms of σx′ and σy′ .

%% c.) Plot and compare the distributions of fxy(x, y) and fx′y′(x′, y′).
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
xlabel('x')
ylabel('y')
zlabel('Probability Density')
title("Positive Correlation f(x,y) ")
xprime_yprime = R_matrix * X.';
mu = 0;
y_decorr = mvnpdf(xprime_yprime',mu,Sigma);
y_decorr = reshape(y_decorr,length(x2),length(x1));
figure(2)
surf(x1,x2,y_decorr)
axis([-3 3 -3 3 0 0.4])
xlabel('x prime')
ylabel('y prime')
zlabel('Probability Density')
title("Uncorrelated f(x,y)")


% Constants
c = 3e8;  % Speed of light (m/s)
max_dimension_handset = 0.15;  % Maximum linear dimension of the handset (m)

% Estimate decorrelation length
decorrelation_length = max_dimension_handset;

% Calculate carrier frequency for which decorrelation length matches the size of the handset
carrier_frequency = c / decorrelation_length;

msg = sprintf("Estimated Carrier Frequency: %s Hz", carrier_frequency);
disp(msg)


gamma = linspace(1e-4,1e5,1e5);
% Define Eb/N0 values in dB
EbN0dB = 0:20; % Range of Eb/N0 values in dB
% Convert Eb/N0 to linear scale
EbN0 = 10.^((EbN0dB / 10));

