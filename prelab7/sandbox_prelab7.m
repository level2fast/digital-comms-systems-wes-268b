% % Generate an example input signal (you can replace this with your actual input signal)
% clear
% clf
% % Parameters
% num_taps = 32;        % Number of taps
% max_delay = 16;       % Maximum delay (samples)
% max_attenuation = .9;  % Maximum attenuation
% 
% % Generate random taps
% taps = zeros(num_taps, 2);
% for i = 1:num_taps
%     %taps(i, 1) = randi(max_delay);           % Random delay
%     taps(i, 1) = max_delay;          % Random delay
%     taps(i, 2) = max_attenuation;    % Random attenuation
% end
% 
% % Create channel impulse response
% channel_length = max(taps(:, 1));
% channel_impulse_response = zeros(1, channel_length);
% channel_impulse_response(1) = 1;
% for i = 1:num_taps
%     delay = taps(i, 1);
%     attenuation = taps(i, 2);
%     channel_impulse_response(delay) = attenuation; % Shift by one to account for MATLAB indexing
% end
% 
% % Normalize channel impulse response to have unit energy
% channel_impulse_response = channel_impulse_response / norm(channel_impulse_response);
% 
% % Display channel impulse response
% disp('Channel Impulse Response:');
% disp(channel_impulse_response);
% 
% % Generate transmitted signal (random symbols)
% num_symbols = 32;
% transmitted_signal = randn(1, num_symbols);
% 
% % Convolve transmitted signal with channel impulse response
% received_signal = conv(transmitted_signal, channel_impulse_response);
% 
% % Display results
% figure(4);
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
% 
% 


% input_signal = randn(1, 32);
% 
% % Apply the channel impulse response
% output_signal = channel_impulse_response(input_signal);
% 
% % Plot the input and output signals
% figure;
% subplot(2,1,1);
% plot(input_signal);
% title('Input Signal');
% subplot(2,1,2);
% plot(output_signal);
% title('Output Signal with Channel Impulse Response');
% 




% % Parameters
% num_taps = 5; % Number of taps
% max_delay = 10; % Maximum delay (in samples)
% max_attenuation = 1; % Maximum attenuation
% 
% % Generate random tap delays and attenuations
% tap_delays = randi(max_delay, 1, num_taps);
% tap_attenuations = rand(1, num_taps) * max_attenuation;
% 
% % Sort tap delays in ascending order
% [tap_delays, sort_idx] = sort(tap_delays);
% tap_attenuations = tap_attenuations(sort_idx);
% 
% % Create tapped delay line impulse response
% h = zeros(1, max_delay + 1); % Initialize impulse response
% for i = 1:num_taps
%     h(tap_delays(i) + 1) = tap_attenuations(i);
% end
% 
% % Plot impulse response
% figure(1)
% stem(0:max_delay,h);
% xlabel('Delay (samples)');
% ylabel('Attenuation');
% title('Tapped Delay Line Impulse Response');
% 
% % Display tap delays and attenuations
% disp('Tap Delays:');
% disp(tap_delays);
% disp('Tap Attenuations:');
% disp(tap_attenuations);
% 
% 
% 
% syms x
% fplot(dirac(0:10))

% channel_mat_h = [complex(-0.82653,0.4552209) complex(-0.413536, +1.0669);
%                  complex(-0.413536,-1.0669) complex(0.82653,0.4552209)];
% [~,eig_vals] = eig(channel_mat_h);
% disp(abs((eig_vals)))


channel_mat_h = [complex(-0.82653,0.4552209) complex(-0.82653,0.4552209);
                 conj(complex(-0.82653,0.4552209)) -conj(complex(-0.82653,0.4552209))];
disp(det(channel_mat_h))
result = channel_mat_h * inv(channel_mat_h);
disp(result)

