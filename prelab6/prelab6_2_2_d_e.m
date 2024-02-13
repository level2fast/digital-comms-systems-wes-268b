%% Prelab6 Problem 2.2.d
% Find an input vector x such that the output vector y expressing the time 
% series is a complex sinusoid at frequency ω = 2pi/N, where N = 8.
[dft_twiddle_factors_mat] = calc_dft_twiddle_factors([1 1 1 1 1 1 1 1]);
[idft_twiddle_factors_mat] = calc_idft_twiddle_factors([1 1 1 1 1 1 1 1]);
omega = 2*pi/8;
y = exp(1j*omega.*(0:7));

figure(1)
subplot(2, 1, 1);
plot(real(y))
title('Real Part of y');
xlabel('time');
ylabel('Amplitude');

disp(y)
disp(fft(y))
x = y * dft_twiddle_factors_mat;
subplot(2, 1, 2);
plot(real(x)/norm(x))
title('Real Part of x');
xlabel('time');
ylabel('Amplitude');


%% Problem 2.2.e
omega = 4*pi/8;
y = cos(omega.*(0:7));
figure(2)
subplot(2, 1, 1);
plot((y))
title('Real Part of y');
xlabel('time');
ylabel('Amplitude');

disp(y)
disp(fft(y))
x = y * dft_twiddle_factors_mat;
subplot(2, 1, 2);
plot(real(x)/norm(x))
title('Real Part of x');
xlabel('time');
ylabel('Amplitude');






