clc; clear; close all;

%Prelab 6 Matlab Simulation

% a) generate a single, QPSK modulated subcarrier, w/ random data using 
%   rectangular baseband pulses. 32 samples per symbol.

symbols = 1024;
samplesPerSymbol = 32;
totalSamples = symbols * samplesPerSymbol;

%Determine Time Range
samplingFrequency = 1e3;
samplingPeriod = 1/samplingFrequency;
totalTimeDuration = (totalSamples-1) * samplingPeriod;
timeRange = 0:samplingPeriod:totalTimeDuration; 

%Determine Frequency Range
nyquistFrequency = samplingFrequency / 2;
frequencySpacing = samplingFrequency/totalSamples;   
frequencyRange = -nyquistFrequency:frequencySpacing:...
                 (nyquistFrequency-frequencySpacing);

%Create a data source which consists of a random sequence of 0s and 1s.
phaseComponentsPerSymbol = 2;
randomSequenceSize = symbols * phaseComponentsPerSymbol;
psuedoNoiseOrder = 5;
randomSequence = pngen(psuedoNoiseOrder, randomSequenceSize);

% BPSK maps 1's to 1 and 0's to -1
bpskRandomSequence = changem(randomSequence, -1, 0);

%Separate the data source into two data streams for the I and Q channels 
inPhaseChannelSequence = bpskRandomSequence(1:2:end);
quadraturePhaseChannelSequence = bpskRandomSequence(2:2:end);

%Recombine the In-Phase and Quad-Phase components to create complex symbols
messageSignal = inPhaseChannelSequence + 1j*quadraturePhaseChannelSequence;

%Create a rectangular baseband pulse
rows = 1;
rectangularPulse = ones(rows,samplesPerSymbol);

%Apply a rectangular baseband pulse to the modulated signal
filteredSignal = kron(messageSignal, rectangularPulse);

%Generate a signal subcarrier signal
getOfdmSubcarrierSignal = @(frequency, time) ...
                            real(exp(-1i*2*pi*frequency.*time));
frequencyDeviation = 0;
subcarrierSignal = getOfdmSubcarrierSignal(frequencyDeviation,timeRange);

%Modulate the message signal with the subcarrier signal
modulatedSignal = filteredSignal.*subcarrierSignal;

%Plot and comment on the power spectrum and time series 
figure('name','Matlab Simulation 1a Time Series of 1 Subcarrier');

subplot(2,1,1);
plot(timeRange,real(modulatedSignal));
title("Time Series of QPSK In-Phase Component");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

subplot(2,1,2);
plot(timeRange,imag(modulatedSignal));
title("Time Series of QPSK Quad-Phase Component");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

getFrequencySpectrum = @(timeWaveform) fftshift(fft(timeWaveform));

getPowerSpectrum = @(timeWaveform) abs(...
                  getFrequencySpectrum(timeWaveform)).^2;

convertLinearTodB = @(powerSpectrum) 10*log10(powerSpectrum);

figure('name','Matlab Simulation 1a Power Spectrum');
plot(frequencyRange,convertLinearTodB(getPowerSpectrum(modulatedSignal)));
title("Power Spectrum of QPSK");
ylabel("Power (dB)")
xlabel("Frequency (Hz)")
grid on;

% b)Generate 2 additional modulated subcarriers with a positive and 
% negative frequency deviation 
frequencyDeviation = samplingFrequency / samplesPerSymbol;
subcarrier1 = getOfdmSubcarrierSignal(frequencyDeviation, timeRange);
subcarrier2 = getOfdmSubcarrierSignal(-frequencyDeviation, timeRange);

%Create a data source which consists of a random sequence of 0s and 1s.
randomSequence1 = pngen(psuedoNoiseOrder, randomSequenceSize);
randomSequence2 = pngen(psuedoNoiseOrder, randomSequenceSize);

% BPSK maps 1's to 1 and 0's to -1
bpskRandomSequence1 = changem(randomSequence1, -1, 0);
bpskRandomSequence2 = changem(randomSequence2, -1, 0);

%Separate the data source into two data streams for the I and Q channels 
inPhaseChannelSequence1 = bpskRandomSequence1(1:2:end);
quadraturePhaseChannelSequence1 = bpskRandomSequence1(2:2:end);

inPhaseChannelSequence2 = bpskRandomSequence2(1:2:end);
quadraturePhaseChannelSequence2 = bpskRandomSequence2(2:2:end);

%Recombine the In-Phase and Quad-Phase components to create complex symbols
messageSignal1 = inPhaseChannelSequence1 ...
                + 1j*quadraturePhaseChannelSequence1;

messageSignal2 = inPhaseChannelSequence2 ...
                + 1j*quadraturePhaseChannelSequence2;

%Apply a rectangular baseband pulse to the modulated signal
filteredSignal1 = kron(messageSignal1, rectangularPulse);
filteredSignal2 = kron(messageSignal2, rectangularPulse);

%Modulate the message signal with the subcarrier signal
modulatedSignal1 = filteredSignal1.*subcarrier1;
modulatedSignal2 = filteredSignal2.*subcarrier2;

%Combine all 3 subcarriers together to create 1 signal
combinedSignal = modulatedSignal + modulatedSignal1 + modulatedSignal2;

%Plot the Time Series of the 3 Subcarriers
figure('name','Matlab Simulation 1b Time Series of 3 Subcarriers');

subplot(2,1,1);
plot(timeRange,real(combinedSignal));
title("Time Series of In-Phase Component");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

subplot(2,1,2);
plot(timeRange,imag(combinedSignal));
title("Time Series of Quad-Phase Component");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

figure('name','Matlab Simulation 1b Power Spectrum');
plot(frequencyRange,convertLinearTodB(getPowerSpectrum(combinedSignal)));
title("Power Spectrum of QPSK");
ylabel("Power (dB)")
xlabel("Frequency (Hz)")
ylim([-10 90]);
xlim([-500 500]);
grid on;

% 1ci)
numSubcarriers = 32;
newRandSeqSize = numSubcarriers * randomSequenceSize;
randomSequence3 = pngen(psuedoNoiseOrder, newRandSeqSize);

bpskRandSeq3 = changem(randomSequence3, -1, 0);

inPhase3 = bpskRandSeq3(1:2:end);
quadPhase3 = bpskRandSeq3(2:2:end);

qpskSymbols = inPhase3 + 1j*quadPhase3;

parallelizerBlocks = zeros(numSubcarriers, symbols);

qpskSignalIndex = 1;

% starting from the subcarrier index, save every 32nd element into a 
% parallel block
for subcarrierIndex = 1:numSubcarriers    
    parallelizerBlocks(subcarrierIndex,:) = qpskSymbols(subcarrierIndex:...
                                            numSubcarriers:end);
end

%Create an IFFT matrix where each element is given by 
% [F^-1]nk = (1/N) exp((1j*2*pi*n*k)/N)
% n = row index
% k = column index
% n & k are 0-based indexed
normalizationFactor = 1/sqrt(numSubcarriers);
ifftFrequency = (2*pi) / samplesPerSymbol;
getIfftElement = @(row,col) (1/samplesPerSymbol) * ...
                            exp(1j*ifftFrequency*row*col);

ifft = zeros(samplesPerSymbol, samplesPerSymbol);

for rowIndex = 1:samplesPerSymbol
    for columnIndex = 1:samplesPerSymbol
        ifft(rowIndex,columnIndex) = ...
            normalizationFactor * getIfftElement(rowIndex-1,columnIndex-1);
    end
end

timeParallel = ifft*parallelizerBlocks;

serializedTime = reshape(timeParallel, 1, []);

figure('name','Matlab Simulation 1di Power Spectrum of 32 Subcarriers');
plot(frequencyRange,convertLinearTodB(getPowerSpectrum(serializedTime)));
title("Power Spectrum of 32 Subcarriers");
ylabel("Power (dB)")
xlabel("Frequency (Hz)")
grid on;

% 1cii
numSubcarriersInSubset = 3;
parallelizerBlockFor3Subcarriers = parallelizerBlocks;
for subcarrierIndex = numSubcarriersInSubset+1:numSubcarriers
    parallelizerBlockFor3Subcarriers(subcarrierIndex,:) = zeros(1,symbols);
end
timeParallel = ifft*parallelizerBlockFor3Subcarriers ;
serializedTimeFor3Subcarriers = reshape(timeParallel, 1, []);

figure('name','Matlab Simulation 1dii Power Spectrum of 32 Subcarriers');
plot(frequencyRange,convertLinearTodB(getPowerSpectrum(...
                    serializedTimeFor3Subcarriers)));
title("Power Spectrum of 3 Subcarriers");
ylabel("Power (dB)")
xlabel("Frequency (Hz)")
grid on;

% 1ciii) When all the subcarrier constellations are overlaid ontop of one
% another, the constellation plot appears more dense packed around (0,0)
% (i.e. In-Phase and Quad-Phase components are at 0)

%If I tried to generate an eye-pattern to find the optimal sampling time, I
%would have to be very precise, otherwise I would fail to properly sample
%the waveform. In fact, I would have to sample in the frequency domain
%instead
figure('name','Matlab Simulation 1ciii Time Series of 32 Subcarriers');
plot(timeRange,real(serializedTime));
title("Time Series of 32 Subcarriers");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

figure('name','Matlab Simulation 1ciii Time Series of 3 Subcarriers');
plot(timeRange,real(serializedTimeFor3Subcarriers));
title("Time Series of 3 Subcarriers");
ylabel("Amplitude")
xlabel("Time (seconds)")
grid on;

figure('name','Matlab Simulation 1ciii Constellation of 32 Subcarriers');
scatter(real(serializedTime), imag(serializedTime))
title("Constellation of 32 Subcarriers");
ylabel("Imaginary")
xlabel("Real")

figure('name','Matlab Simulation 1ciii Constellation of 3 Subcarriers');
scatter(real(serializedTimeFor3Subcarriers), ...
    imag(serializedTimeFor3Subcarriers))
title("Constellation of 3 Subcarriers");
ylabel("Imaginary")
xlabel("Real")

% d Noise
mean = 0;
standardDeviation = 0.05;
inPhaseNoise = normrnd(mean, standardDeviation, rows, totalSamples);
quadPhaseNoise = normrnd(mean, standardDeviation, rows, totalSamples);
noise = inPhaseNoise + 1j*quadPhaseNoise;
noisySignal = serializedTime + noise;

% e Noisy Recieved Values using FFT
%The fft is the complex conjugate of the ifft
getfftElement = @(row,col) (1/samplesPerSymbol) * ...
                            exp(-1j*ifftFrequency*row*col);

fft = zeros(samplesPerSymbol, samplesPerSymbol);

for rowIndex = 1:samplesPerSymbol
    for columnIndex = 1:samplesPerSymbol
        fft(rowIndex,columnIndex) = ...
            normalizationFactor * getfftElement(rowIndex-1,columnIndex-1);
    end
end

recievedNoisySignal = fft*reshape(noisySignal, samplesPerSymbol, []);

%f QPSK Constellations based on frequency
figure('name',['Matlab Simulation 1f Constellation of Highest ...' ...
    'Subcarrier Frequency']);
scatter(real(recievedNoisySignal(numSubcarriers/2,:)), ...
    imag(recievedNoisySignal(numSubcarriers/2,:))         )
title("Constellation of the Highest Frequency Subcarrier");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1f Constellation of the Lowest ...' ...
    'Subcarrier Frequency']);
scatter(real(recievedNoisySignal(2,:)), ...
    imag(recievedNoisySignal(2,:)))
title("Constellation of the Lowest Frequency Subcarrier");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1f Constellation of the ...' ...
    'Subcarrier Frequency Halfway b/t the Highest and Lowest']);
scatter(real(recievedNoisySignal(3*numSubcarriers/4,:)), ...
    imag(recievedNoisySignal(3*numSubcarriers/4,:)))
title("Constellation of the Subcarrier Frequency Halfway b/t the ..." + ...
    "Highest and Lowest");
ylabel("Imaginary")
xlabel("Real")

%g
offsetSignal = [serializedTime(3:end) serializedTime(1:2)];
offsetSignal = reshape(offsetSignal, numSubcarriers, []);

recievedOffset = fft*offsetSignal;

figure('name',['Matlab Simulation 1g']);
scatter(real(recievedOffset(numSubcarriers/2,:)), ...
    imag(recievedOffset(numSubcarriers/2,:))         )
title("Constellation of the Highest Frequency Subcarrier w/ a 2 ..." + ...
    "sample offset");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1g']);
scatter(real(recievedOffset(2,:)), ...
    imag(recievedOffset(2,:)))
title("Constellation of the Lowest Frequency Subcarrier w/ a 2 ..." + ...
    "Sample Offset");
ylabel("Imaginary")
xlabel("Real")

figure('name',['Matlab Simulation 1g']);
scatter(real(recievedOffset(3*numSubcarriers/4,:)), ...
    imag(recievedOffset(3*numSubcarriers/4,:)))
title("Constellation of the Middle Subcarrier w/ a 2 sample offset");
ylabel("Imaginary")
xlabel("Real")

% gi) With every symbol offset, the subcarrier constellations of the
% highest, lowest, and middle frequency subcarriers appear to shift in
% phase. 