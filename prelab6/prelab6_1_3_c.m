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
% disp(['Integral value: ', num2str(y_t)]);
% figure(1)
% plot(delta_f_vec,y_t)
% title(['Symbol Energy of $\cos({2} \pi f_c \tau) \lbrack 2\cos({2} \pi ' ...
%     '(f_c + \Delta f) \tau) \rbrack in \Delta f$ steps'],'Interpreter','latex');
% xlabel('$\Delta f$','Interpreter','latex');
% ylabel('Symbol Energy y(t)');
% grid on;

%% Plot BER as a function of delta_f
% Need to create a new vector has values for BER based on the delta_f value
% set Eb/n0=5 should be deltaf= 0 ,
% This is a factor of 2 away from -3dB
%Initial value for delta_f zero should be around 10-e2, sweep deltaf to the right to see what the plot looks like
Start_Eb_N0_dB = db2pow(5);
Eb_N0_dB = Start_Eb_N0_dB:1:14;                % Range of Eb/N0 values in dB
Eb_N0    = 10.^(Eb_N0_dB / 10);    % Convert dB to linear scale
BER      = 0.5* qfunc(Eb_N0);
delta_f_sweep = 0:10e-2:10e-1;

% Calculate the theoretical Pe for BPSK
% Pe_theoretical = qfunc(sqrt(Eb_N0));

% Find the formula for theoritical Pe as delta f varies
% Plot the results
figure(1);
semilogy(Eb_N0, delta_f_sweep, 'LineWidth', 2);
title('BER vs. delta_f for OFDM');
xlabel('$ \Delta f (Initial \Delta f = 0)$','Interpreter','latex');
ylabel('BER (Initial EbN0 = 5 linear)');
grid on;