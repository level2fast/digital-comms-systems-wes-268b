%% Generation of OFDM Signals

clear
% Parameters
num_subcarriers = 4;              % Number of subcarriers
f_0 = 1;                          % center frequency 
bw_hz = 2;                        % bandwidth
delta_f  = bw_hz/num_subcarriers; % subcarrier spacing
Ts_sec = 1/delta_f;               % symbol duration
symbol_rate_sps = 1/(Ts_sec*50);  % symbol rate, i.e. the rate at which symbols are transmitted within each subcarrier
t = 1:symbol_rate_sps:Ts_sec;     % time vector

% Generate complex exponential carriers
n_subcarrier = (1:num_subcarriers)';
rows = size(n_subcarrier,1);
cols = size(t,2);
s_t = zeros(rows,cols);

%% a) Plot the time waveform for each of the four subcarriers using f0 = 1.
for nIdx = 1:size(n_subcarrier)
    sub_carrier = n_subcarrier(nIdx);
    s_t(nIdx,:) = exp(-1j * 2 * pi * sub_carrier * f_0 * t);
    % Plot the generated OFDM signal
    figure(1);
    subplot(2, 2, nIdx);
    plot(t,(s_t(nIdx,:)));
    title('Real Part of OFDM Signal');
    xlabel('time');
    ylabel('Amplitude');
    legend("OFDM real part")
    grid on
    hold on
end
hold off
%% b) Plot the spectrum for each of the four subcarriers on the same plot.
for nIdx = 1:size(n_subcarrier)
    % Plot the generated OFDM signal
    figure(2);
    plot(t,real(s_t(nIdx,:)));
    title('Real Part of OFDM Signal');
    xlabel('time');
    ylabel('Amplitude');
    legend
    grid on
    hold on
end
hold off



