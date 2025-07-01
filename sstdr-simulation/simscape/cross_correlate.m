
% Take in signal
time = out.simout.Time;
measurement = out.simout.Data(:);
code = code_interp;

flip_code = flipud(code);
code_padded = [code; zeros(length(measurement) - length(code), 1)];

corr1 = conv(measurement, flip_code);
corr2 = xcorr(measurement, code);

close all;

figure(1);
plot(corr1');

% figure(4);
% plot(corr2');

figure(2);
plot(time, code_padded, 'red');
hold on;
plot(time, measurement, 'blue');

% Plot difference between measurement and code
figure(5);
difference = measurement - code_padded;
plot(difference, 'magenta');
title('Difference (Measurement - Code) - Code Zero-Padded');

xlabel('Sample');
ylabel('Amplitude');
grid on;
