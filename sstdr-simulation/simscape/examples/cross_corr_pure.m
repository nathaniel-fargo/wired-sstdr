%% SSTDR Cross-Correlation Simulation with Modulated PN Sequence
% This script demonstrates spread spectrum time domain reflectometry (SSTDR)
% by generating a modulated PN sequence, simulating a reflection, and 
% performing cross-correlation analysis

clear; close all; clc;

%% Generate PN Sequence using gen_pn_code.m
fprintf('Generating PN sequence using gen_pn_code.m...\n');
run('gen_pn_code.m');

% Get the generated PN sequence from base workspace
pn_seq = evalin('base', 'code_seq');
pn_seq = pn_seq(:)'; % Ensure it's a row vector
pn_length = length(pn_seq);

fprintf('Using %d-bit PN sequence with length %d\n', log2(pn_length+1), pn_length);

%% Parameters based on PN sequence
fs = 1e6;               % Sampling frequency (Hz) - same as gen_pn_code.m
T = pn_length / fs;     % Total time duration based on PN length
t = (0:pn_length-1) / fs; % Time vector based on PN sequence (row vector)
fc = 100e3;             % Carrier frequency for modulation (Hz)
reflection_delay_samples = round(0.3 * pn_length); % Reflection delay in samples
reflection_amplitude = -0.5; % Amplitude of reflection (negative for inversion)

fprintf('Sequence duration: %.3f ms\n', T*1000);
fprintf('Reflection delay: %d samples (%.3f ms)\n', reflection_delay_samples, reflection_delay_samples/fs*1000);

%% Modulate PN sequence with sine wave
carrier = sin(2*pi*fc*t);
modulated_pn = pn_seq .* carrier;

%% Create measurement trace
measurement = modulated_pn; % Direct signal at time 0

% Add reflection at specified delay
if reflection_delay_samples < length(t)
    reflection_signal = zeros(size(t));
    end_idx = min(length(modulated_pn), length(t) - reflection_delay_samples);
    reflection_signal((reflection_delay_samples+1):(reflection_delay_samples+end_idx)) = ...
        reflection_amplitude * modulated_pn(1:end_idx);
    measurement = measurement + reflection_signal;
end

%% Add some noise to make it more realistic
noise_level = 0.05;
measurement = measurement + noise_level * randn(size(measurement));

%% Perform cross-correlation
fprintf('Performing cross-correlation...\n');
% Ensure both signals are row vectors for xcorr
measurement = measurement(:)';
modulated_pn = modulated_pn(:)';
[correlation, lags] = xcorr(measurement, modulated_pn);
time_lags = lags / fs;

% Normalize correlation
correlation = correlation / max(abs(correlation));

%% Frequency Domain Cross-Correlation
fprintf('Performing frequency domain cross-correlation...\n');

% Zero-pad both signals to avoid circular correlation artifacts
N_fft = 2^nextpow2(2*length(measurement) - 1);

% Take FFTs
X = fft(measurement, N_fft);
Y = fft(modulated_pn, N_fft);

% Cross-correlation in frequency domain: IFFT(X * conj(Y))
correlation_freq = ifft(X .* conj(Y));

% The FFT-based correlation needs proper arrangement to match xcorr
% After IFFT, we have a circular correlation that needs to be converted to linear
N = length(modulated_pn);

% Extract the linear correlation portion (first 2*N-1 points)
correlation_freq = correlation_freq(1:(2*N-1));

% Rearrange to match xcorr output: move second half to beginning
% xcorr arranges as [negative lags, zero lag, positive lags]
% FFT result: [lag 0, lag 1, ..., lag N-1, lag -N+1, ..., lag -1]
% We want: [lag -N+1, ..., lag -1, lag 0, lag 1, ..., lag N-1]
correlation_freq = [correlation_freq(N+1:end), correlation_freq(1:N)];

% Create matching lag vector
max_lag = N - 1;
lags_freq = -max_lag:max_lag;
time_lags_freq = lags_freq / fs;

% Normalize
correlation_freq = correlation_freq / max(abs(correlation_freq));

% Trim to same range as time domain correlation for comparison
valid_range = (time_lags_freq >= min(time_lags)) & (time_lags_freq <= max(time_lags));
correlation_freq_trimmed = correlation_freq(valid_range);
time_lags_freq_trimmed = time_lags_freq(valid_range);

%% Find peaks in correlation
[peaks, peak_locs] = findpeaks(abs(correlation), 'MinPeakHeight', 0.1, 'MinPeakDistance', 50);
peak_times = time_lags(peak_locs);

fprintf('Detected peaks at times: ');
fprintf('%.3f ', peak_times);
fprintf('seconds\n');

%% Plot results
figure('Position', [100, 100, 1400, 1000]);

% Original PN sequence
subplot(3,2,1);
plot(t, pn_seq, 'b-', 'LineWidth', 1);
title('Original PN Sequence (Before Modulation)', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0, min(0.1e-3, T/10)]); % Show first part for clarity

% Modulated PN sequence
subplot(3,2,2);
plot(t, modulated_pn, 'r-', 'LineWidth', 1);
title('Modulated PN Sequence (Carrier Modulated)', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0, min(0.1e-3, T/10)]); % Show first part for clarity

% Measurement trace
subplot(3,2,3);
plot(t, measurement, 'g-', 'LineWidth', 1);
title('Measurement Trace (Direct + Reflection + Noise)', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Time domain cross-correlation result
subplot(3,2,4);
plot(time_lags, correlation, 'k-', 'LineWidth', 1.5);
hold on;
% Mark detected peaks
if ~isempty(peak_times)
    plot(peak_times, correlation(peak_locs), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    % Add text annotations for peaks
    for i = 1:length(peak_times)
        text(peak_times(i), correlation(peak_locs(i)) + 0.1, ...
             sprintf('%.3fs', peak_times(i)), 'HorizontalAlignment', 'center');
    end
end
title('Time Domain Cross-Correlation', 'FontSize', 12);
xlabel('Time Lag (s)');
ylabel('Normalized Correlation');
grid on;
xlim([-T/2, T*1.5]);
ylim([-0.8, 1.1]);

% Add vertical lines at expected locations
line([0, 0], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);
reflection_delay_time = reflection_delay_samples / fs;
line([reflection_delay_time, reflection_delay_time], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
legend('Correlation', 'Detected Peaks', 'Direct Signal', 'Expected Reflection', 'Location', 'best');

% Frequency domain cross-correlation result
subplot(3,2,5);
plot(time_lags_freq_trimmed, real(correlation_freq_trimmed), 'b-', 'LineWidth', 1.5);
hold on;
% Find peaks in frequency domain correlation too
[peaks_freq, peak_locs_freq] = findpeaks(abs(correlation_freq_trimmed), 'MinPeakHeight', 0.1, 'MinPeakDistance', 50);
if ~isempty(peak_locs_freq)
    peak_times_freq = time_lags_freq_trimmed(peak_locs_freq);
    plot(peak_times_freq, real(correlation_freq_trimmed(peak_locs_freq)), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
end
title('Frequency Domain Cross-Correlation', 'FontSize', 12);
xlabel('Time Lag (s)');
ylabel('Normalized Correlation');
grid on;
xlim([-T/2, T*1.5]);
ylim([-0.8, 1.1]);
% Add vertical lines at expected locations
line([0, 0], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);
line([reflection_delay_time, reflection_delay_time], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
legend('Correlation', 'Detected Peaks', 'Direct Signal', 'Expected Reflection', 'Location', 'best');

% Comparison of time vs frequency domain correlation
subplot(3,2,6);
plot(time_lags, correlation, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Time Domain');
hold on;
plot(time_lags_freq_trimmed, real(correlation_freq_trimmed), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Frequency Domain');
title('Time vs Frequency Domain Comparison', 'FontSize', 12);
xlabel('Time Lag (s)');
ylabel('Normalized Correlation');
grid on;
xlim([-T/2, T*1.5]);
ylim([-0.8, 1.1]);
% Add vertical lines at expected locations
line([0, 0], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);
line([reflection_delay_time, reflection_delay_time], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
legend('show', 'Location', 'best');

%% Display results
fprintf('\n=== SSTDR Analysis Results ===\n');
fprintf('Sampling frequency: %.0f Hz\n', fs);
fprintf('Carrier frequency: %.0f Hz\n', fc);
fprintf('PN sequence length: %d chips\n', pn_length);
fprintf('Expected reflection delay: %.3f us\n', reflection_delay_time*1e6);
fprintf('Expected reflection amplitude: %.2f\n', reflection_amplitude);
fprintf('Noise level: %.3f\n', noise_level);

if length(peak_times) >= 2
    measured_delay = abs(peak_times(2) - peak_times(1));
    fprintf('Measured reflection delay: %.3f us\n', measured_delay*1e6);
    delay_error = abs(measured_delay - reflection_delay_time);
    fprintf('Delay measurement error: %.3f us (%.1f%%)\n', delay_error*1e6, 100*delay_error/reflection_delay_time);
end

fprintf('\nThe plot shows:\n');
fprintf('1. Original PN sequence (binary)\n');
fprintf('2. Sine wave modulated PN sequence\n');
fprintf('3. Simulated measurement with direct signal + reflection + noise\n');
fprintf('4. Time domain cross-correlation result showing sinc-like peaks\n');
fprintf('5. Frequency domain cross-correlation result (should match time domain)\n');
fprintf('6. Direct comparison between time and frequency domain methods\n');
fprintf('   - Peak at t=0: Direct signal correlation\n');
fprintf('   - Peak at t=%.3fus: Reflection correlation (inverted)\n', reflection_delay_time*1e6);

%% Compare computational efficiency
fprintf('\n=== Method Comparison ===\n');
fprintf('Time domain correlation: Uses built-in xcorr() function\n');
fprintf('Frequency domain correlation: Uses FFT-based approach\n');
fprintf('- Zero-padded to %d points for linear correlation\n', N_fft);
fprintf('- Frequency domain is typically faster for long sequences\n');
fprintf('- Both methods should produce nearly identical results\n');

% Calculate correlation between the two methods as a validation
if length(correlation) == length(correlation_freq_trimmed)
    method_correlation = corrcoef(correlation, real(correlation_freq_trimmed));
    fprintf('Correlation between methods: %.6f (should be ~1.0)\n', method_correlation(1,2));
else
    fprintf('Note: Different vector lengths prevent direct method comparison\n');
end
