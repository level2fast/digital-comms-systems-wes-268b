function output_signal = convolve_with_channel(cir)
    % Define impulse response with a delay of N samples
    arguments
        cir.input_signal (1,:) = 0
        cir.time_delay (1,1) = 16
    end
    input_signal = cir.input_signal;
    time_delay   = length(input_signal)/2;
    
    % first create the channel impulse response
    impulse_response = [zeros(1, time_delay), 1];

    % Convolve input signal with impulse response
    output_signal = conv(input_signal, impulse_response);

    % Truncate the output signal to match the length of the input signal
    output_signal = output_signal(1:length(input_signal));
end