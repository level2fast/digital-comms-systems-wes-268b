% Parameters
num_taps = 5; % Number of taps
max_delay = 10; % Maximum delay (in samples)
max_attenuation = 1; % Maximum attenuation

% Generate random tap delays and attenuations
tap_delays = randi(max_delay, 1, num_taps);
tap_attenuations = rand(1, num_taps) * max_attenuation;

% Sort tap delays in ascending order
[tap_delays, sort_idx] = sort(tap_delays);
tap_attenuations = tap_attenuations(sort_idx);

% Create tapped delay line impulse response
h = zeros(1, max_delay + 1); % Initialize impulse response
for i = 1:num_taps
    h(tap_delays(i) + 1) = tap_attenuations(i);
end

% Plot impulse response
figure(1)
stem(0:max_delay,h);
xlabel('Delay (samples)');
ylabel('Attenuation');
title('Tapped Delay Line Impulse Response');

% Display tap delays and attenuations
disp('Tap Delays:');
disp(tap_delays);
disp('Tap Attenuations:');
disp(tap_attenuations);



syms x
fplot(dirac(0:10))
