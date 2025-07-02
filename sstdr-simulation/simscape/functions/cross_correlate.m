function [results] = cross_correlate(varargin)
%CROSS_CORRELATE Enhanced cross-correlation analysis for SSTDR simulations
%
% Usage:
%   cross_correlate()                               % Use simulation output from base workspace
%   cross_correlate('method', 'freq')               % Use frequency domain
%   cross_correlate('method', 'time')               % Use time domain
%   cross_correlate('method', 'both')               % Compare both methods
%   results = cross_correlate(...)                  % Return results structure
%
% Parameters:
%   'method'        - 'time', 'freq', 'both' (default: 'freq')
%   'sim_output'    - Simulation output structure (default: 'out' from base)  
%   'reference_var' - Reference code variable name (default: 'pn_interp_modulated')
%   'plot_results'  - Show plots (default: true)
%   'find_peaks'    - Enable peak detection (default: true)
%   'peak_threshold'- Minimum peak height (default: 0.1)
%   'normalize'     - Normalize correlation (default: true)
%   'export_results'- Export to base workspace (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'method', 'freq', @(x) ismember(x, {'time', 'freq', 'both'}));
addParameter(p, 'sim_output', [], @isstruct);
addParameter(p, 'reference_var', 'pn_interp_modulated', @ischar);
addParameter(p, 'plot_results', true, @islogical);
addParameter(p, 'find_peaks', true, @islogical);
addParameter(p, 'peak_threshold', 0.1, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'normalize', true, @islogical);
addParameter(p, 'export_results', true, @islogical);

parse(p, varargin{:});
cfg = p.Results;

%% Get simulation data
if isempty(cfg.sim_output)
    try
        sim_data = evalin('base', 'out');
    catch
        error('No simulation output found. Please provide sim_output or ensure "out" exists in base workspace.');
    end
else
    sim_data = cfg.sim_output;
end

% Extract measurement signal
time = sim_data.simout.time;
measurement = sim_data.simout.signals.values;
fs = 1 / (time(2) - time(1));  % Calculate sampling frequency

fprintf('Processing simulation data:\n');
fprintf('  Duration: %.3f ms\n', max(time)*1000);
fprintf('  Samples: %d\n', length(measurement));
fprintf('  Sampling rate: %.1f kHz\n', fs/1000);

%% Get reference code
try
    if evalin('base', sprintf('exist(''%s'', ''var'')', cfg.reference_var))
        reference_code = evalin('base', cfg.reference_var);
        fprintf('Using reference code: %s (%d samples)\n', cfg.reference_var, length(reference_code));
    else
        % Fallback to legacy variable
        reference_code = evalin('base', 'code_interp');
        fprintf('Using fallback reference code: code_interp (%d samples)\n', length(reference_code));
    end
catch
    error('Reference code not found. Please run gen_pn_code() first.');
end

% Ensure column vectors
measurement = measurement(:);
reference_code = reference_code(:);

%% Perform correlation analysis
results = struct();
results.config = cfg;
results.time = time;
results.measurement = measurement;
results.reference_code = reference_code;
results.fs = fs;

if strcmp(cfg.method, 'time') || strcmp(cfg.method, 'both')
    fprintf('Performing time domain cross-correlation...\n');
    tic;
    [correlation_time, lags_time] = xcorr(measurement, reference_code);
    time_lags_time = lags_time / fs;
    compute_time_time = toc;
    
    if cfg.normalize
        correlation_time = correlation_time / max(abs(correlation_time));
    end
    
    results.correlation_time = correlation_time;
    results.time_lags_time = time_lags_time;
    results.compute_time_time = compute_time_time;
    
    fprintf('  Time domain completed in %.3f ms\n', compute_time_time*1000);
end

if strcmp(cfg.method, 'freq') || strcmp(cfg.method, 'both')
    fprintf('Performing frequency domain cross-correlation...\n');
    tic;
    
    % Zero-pad both signals to avoid circular correlation artifacts
    N = max(length(measurement), length(reference_code));
    N_fft = 2^nextpow2(2*N - 1);
    
    % Take FFTs
    X = fft(measurement, N_fft);
    Y = fft(reference_code, N_fft);
    
    % Cross-correlation in frequency domain: IFFT(X * conj(Y))
    correlation_freq = ifft(X .* conj(Y));
    
    % Extract the linear correlation portion and arrange like xcorr
    N_ref = length(reference_code);
    correlation_freq = correlation_freq(1:(2*N_ref-1));
    correlation_freq = [correlation_freq(N_ref+1:end); correlation_freq(1:N_ref)];
    
    % Create time lag vector
    max_lag = N_ref - 1;
    lags_freq = (-max_lag:max_lag)';
    time_lags_freq = lags_freq / fs;
    
    compute_time_freq = toc;
    
    if cfg.normalize
        correlation_freq = correlation_freq / max(abs(correlation_freq));
    end
    
    results.correlation_freq = correlation_freq;
    results.time_lags_freq = time_lags_freq;
    results.compute_time_freq = compute_time_freq;
    results.fft_size = N_fft;
    
    fprintf('  Frequency domain completed in %.3f ms (FFT size: %d)\n', ...
            compute_time_freq*1000, N_fft);
end

%% Peak detection
if cfg.find_peaks
    fprintf('Detecting correlation peaks...\n');
    
    % Calculate MinPeakDistance with safety limits
    % MATLAB findpeaks has a limit of 32735 for MinPeakDistance
    min_peak_distance = min(round(fs/10000), 32000);  % Safe limit
    min_peak_distance = max(min_peak_distance, 1);     % At least 1
    
    if strcmp(cfg.method, 'time') || strcmp(cfg.method, 'both')
        try
            [peaks_time, peak_locs_time] = findpeaks(abs(correlation_time), ...
                'MinPeakHeight', cfg.peak_threshold, 'MinPeakDistance', min_peak_distance);
            peak_times_time = time_lags_time(peak_locs_time);
            
            results.peaks_time = peaks_time;
            results.peak_times_time = peak_times_time;
            results.peak_locs_time = peak_locs_time;
            
            fprintf('  Time domain: %d peaks detected\n', length(peaks_time));
        catch ME
            fprintf('  ⚠ Time domain peak detection failed: %s\n', ME.message);
            fprintf('  Continuing without peak detection...\n');
            results.peaks_time = [];
            results.peak_times_time = [];
            results.peak_locs_time = [];
        end
    end
    
    if strcmp(cfg.method, 'freq') || strcmp(cfg.method, 'both')
        try
            [peaks_freq, peak_locs_freq] = findpeaks(abs(correlation_freq), ...
                'MinPeakHeight', cfg.peak_threshold, 'MinPeakDistance', min_peak_distance);
            peak_times_freq = time_lags_freq(peak_locs_freq);
            
            results.peaks_freq = peaks_freq;
            results.peak_times_freq = peak_times_freq;
            results.peak_locs_freq = peak_locs_freq;
            
            fprintf('  Frequency domain: %d peaks detected\n', length(peaks_freq));
        catch ME
            fprintf('  ⚠ Frequency domain peak detection failed: %s\n', ME.message);
            fprintf('  Continuing without peak detection...\n');
            results.peaks_freq = [];
            results.peak_times_freq = [];
            results.peak_locs_freq = [];
        end
    end
end

%% Performance comparison
if strcmp(cfg.method, 'both')
    speedup = results.compute_time_time / results.compute_time_freq;
    fprintf('\nPerformance comparison:\n');
    fprintf('  Time domain: %.3f ms\n', results.compute_time_time*1000);
    fprintf('  Frequency domain: %.3f ms\n', results.compute_time_freq*1000);
    if speedup > 1
        speedup_text = '(freq faster)';
    else
        speedup_text = '(time faster)';
    end
    fprintf('  Speedup: %.1fx %s\n', abs(speedup), speedup_text);
    
    results.speedup = speedup;
    
    % Verify methods match
    if length(results.correlation_time) == length(results.correlation_freq)
        method_corr = corrcoef(results.correlation_time, real(results.correlation_freq));
        results.method_correlation = method_corr(1,2);
        fprintf('  Method correlation: %.6f (should be ~1.0)\n', results.method_correlation);
    end
end

%% Plotting
if cfg.plot_results
    plot_correlation_results(results, cfg);
end

%% Export results
if cfg.export_results
    fprintf('Exporting results to base workspace...\n');
    assignin('base', 'correlation_results', results);
    
    % Export individual results for backward compatibility
    if isfield(results, 'correlation_time')
        assignin('base', 'corr_time', results.correlation_time);
        assignin('base', 'time_lags', results.time_lags_time);
    end
    if isfield(results, 'correlation_freq')
        assignin('base', 'corr_freq', results.correlation_freq);
        assignin('base', 'freq_lags', results.time_lags_freq);
    end
end

fprintf('Cross-correlation analysis complete!\n');

end

%% Helper function for plotting
function plot_correlation_results(results, cfg)
    
    switch cfg.method
        case 'time'
            figure('Position', [100, 100, 1200, 800]);
            
            subplot(2,1,1);
            plot(results.time, results.measurement, 'b-');
            title('Measurement Signal');
            xlabel('Time (s)');
            ylabel('Amplitude');
            grid on;
            
            subplot(2,1,2);
            plot(results.time_lags_time, results.correlation_time, 'k-', 'LineWidth', 1.5);
            hold on;
            if isfield(results, 'peak_times_time')
                plot(results.peak_times_time, results.correlation_time(results.peak_locs_time), ...
                     'ro', 'MarkerSize', 8, 'LineWidth', 2);
            end
            title('Time Domain Cross-Correlation');
            xlabel('Time Lag (s)');
            ylabel('Correlation');
            grid on;
            
        case 'freq'
            figure('Position', [100, 100, 1200, 800]);
            
            subplot(2,1,1);
            plot(results.time, results.measurement, 'b-');
            title('Measurement Signal');
            xlabel('Time (s)');
            ylabel('Amplitude');
            grid on;
            
            subplot(2,1,2);
            plot(results.time_lags_freq, real(results.correlation_freq), 'r-', 'LineWidth', 1.5);
            hold on;
            if isfield(results, 'peak_times_freq')
                plot(results.peak_times_freq, real(results.correlation_freq(results.peak_locs_freq)), ...
                     'ro', 'MarkerSize', 8, 'LineWidth', 2);
            end
            title('Frequency Domain Cross-Correlation');
            xlabel('Time Lag (s)');
            ylabel('Correlation');
            grid on;
            
        case 'both'
            figure('Position', [100, 100, 1400, 1000]);
            
            subplot(2,2,1);
            plot(results.time, results.measurement, 'b-');
            title('Measurement Signal');
            xlabel('Time (s)');
            ylabel('Amplitude');
            grid on;
            
            subplot(2,2,2);
            ref_time = (0:length(results.reference_code)-1) / results.fs;
            plot(ref_time, results.reference_code, 'g-');
            title('Reference Code');
            xlabel('Time (s)');
            ylabel('Amplitude');
            grid on;
            
            subplot(2,2,3);
            plot(results.time_lags_time, results.correlation_time, 'k-', 'LineWidth', 1.5);
            hold on;
            if isfield(results, 'peak_times_time')
                plot(results.peak_times_time, results.correlation_time(results.peak_locs_time), ...
                     'ro', 'MarkerSize', 6);
            end
            title(sprintf('Time Domain (%.1f ms)', results.compute_time_time*1000));
            xlabel('Time Lag (s)');
            ylabel('Correlation');
            grid on;
            
            subplot(2,2,4);
            plot(results.time_lags_freq, real(results.correlation_freq), 'b-', 'LineWidth', 1.5);
            hold on;
            if isfield(results, 'peak_times_freq')
                plot(results.peak_times_freq, real(results.correlation_freq(results.peak_locs_freq)), ...
                     'ro', 'MarkerSize', 6);
            end
            title(sprintf('Frequency Domain (%.1f ms)', results.compute_time_freq*1000));
            xlabel('Time Lag (s)');
            ylabel('Correlation');
            grid on;
    end
end
