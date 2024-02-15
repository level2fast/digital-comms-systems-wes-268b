%% Prelab6 Problem 3.c For T = 1, determine the maximum value of ∆fmax of the 
%% frequency offset such that the bit error rate does not change by more than 
%% a factor of two when Eb/N0 = 5 (linear).

T = 1; % symbol duration
fc = 1e3;
% Define the function to integrate
func = @(t,Fc,Fd) cos(2*pi*Fc*t) .* (2*cos(2*pi*(Fc + Fd)*t));

delta_f_vec = 0:0.0005:1/(2*T);
y_t = zeros(size(delta_f_vec,2),1)';
for fdIdx = 1:size(delta_f_vec,2)
    % Perform the integration
    y_t(fdIdx) = (1/T) * integral(@(t) func(t,fc,delta_f_vec(fdIdx)), 0, T);
end

% Display the result,
disp(['Integral value: ', num2str(y_t)]);
figure(1)
plot(delta_f_vec,y_t)
title(['BER of $\cos({2} \pi f_c \tau) \lbrack 2\cos({2} \pi ' ...
    '(f_c + \Delta f) \tau$) \rbrack'],'Interpreter','latex');
xlabel('$\Delta f$','Interpreter','latex');
ylabel('BER(Pb)');
grid on;
