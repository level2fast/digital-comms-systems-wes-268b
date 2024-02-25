function [channel_impulse_response] = create_ch_imp_rsp(cir)
%CREATE_CIR Summary of this function goes here
%   Detailed explanation goes here
arguments
    cir.taps       (1,1) = 32
    cir.time_delay (1,1) = 16
end
% Parameters
num_taps = cir.taps;       % Number of taps
max_delay = cir.time_delay; % Maximum delay (samples)
max_attenuation = 0.9;      % Maximum attenuation

% Generate random taps
taps = zeros(num_taps, 2);
for i = 1:num_taps
    taps(i, 1) = (max_delay);     % Random delay
    taps(i, 2) = max_attenuation; % Random attenuation
end

% Create channel impulse response
channel_length = max(taps(:, 1)) + 1;
channel_impulse_response = zeros(1, channel_length);
for i = 1:num_taps
    delay = taps(i, 1);
    attenuation = taps(i, 2);
    channel_impulse_response(delay + 1) = attenuation; % Shift by one to account for MATLAB indexing
end

% Normalize channel impulse response to have unit energy
channel_impulse_response = channel_impulse_response / norm(channel_impulse_response);

% % Display channel impulse response
% disp('Channel Impulse Response:');
% disp(channel_impulse_response);
% 
% % Generate transmitted signal (random symbols)
% num_symbols = 100;
% transmitted_signal = randn(1, num_symbols);
% 
% % Convolve transmitted signal with channel impulse response
%received_signal = conv(transmitted_signal, channel_impulse_response);
%
% Display results
% figure;
% subplot(2, 1, 1);
% stem(0:length(channel_impulse_response)-1, channel_impulse_response, 'filled');
% title('Channel Impulse Response');
% xlabel('Sample');
% ylabel('Amplitude');
%
% subplot(2, 1, 2);
% plot(0:length(received_signal)-1, received_signal);
% hold on;
% stem(0:length(transmitted_signal)-1, transmitted_signal, 'r', 'filled');
% hold off;
% title('Received Signal (with Transmitted Signal)');
% xlabel('Sample');
% ylabel('Amplitude');
% legend('Received Signal', 'Transmitted Signal');
end

