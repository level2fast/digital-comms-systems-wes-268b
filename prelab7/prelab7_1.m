%% 1.a Use the 32 subcarrier data from the previous simulation problem in Prelab 6
% Summary: Delay spread refers to the spreading of a transmitted signal over 
% time due to multipath propagation in the channel. As a result, different 
% copies of the transmitted signal arrive at the receiver with different 
% delays, causing distortion of the received signal. The delay spread can 
% be characterized by the spread of the impulse response of the channel.
%
% Negative effects in OFDM:
% Delay spread can cause intersymbol interference (ISI) and degrade 
% system performance.
%
% How to model:
% One common model used to describe delay spread in OFDM systems is the 
% tapped delay line model. In this model, the channel impulse response is 
% represented as a sum of delayed and attenuated copies of the transmitted 
% signal, referred to as taps. Each tap corresponds to a distinct path in
% the multipath channel.

% Simulate the channel given by the discrete impulse response h[n] = δ[n] − δ[n − 16].
channel_impulse_response = 0;

% Create a delayed version of the received signal. This simulates the
% effect that the channel impulse response has on the transmitted OFDM
% signal.

% Convolution can be used to create a delayed and shifted version of a
% signal. Convolve the transmitted signal with a discrete time impulse
% response.
%% 1.a.i Plot the power density spectrum for the received signal given this channel model



%% 1.a.ii Plot the constellation of several subcarriers 


%% 1.b Now add cyclic prefix(CP) to the data that is longer than the delay in the channel model