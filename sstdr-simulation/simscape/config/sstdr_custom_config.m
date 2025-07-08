function [config] = sstdr_custom_config(varargin)
%SSTDR_CUSTOM_CONFIG Create custom SSTDR configuration with specified frequencies
%
% Usage:
%   sstdr_custom_config()                                    % Interactive mode
%   sstdr_custom_config('carrier_freq', 500e3)               % 500 kHz carrier
%   sstdr_custom_config('chip_rate', 200e3)                  % 200 kHz chip rate
%   sstdr_custom_config('fs', 5e6)                          % 5 MHz sampling
%   sstdr_custom_config('carrier_freq', 1e6, 'chip_rate', 1e6, 'fs', 4e6) % All frequencies
%   config = sstdr_custom_config(...)                       % Return config
%
% Parameters:
%   'carrier_freq'  - Carrier frequency in Hz (default: 100e3)
%   'chip_rate'     - Chip rate in Hz (default: 100e3)
%   'fs'           - Sampling frequency in Hz (default: 400e3)  
%   'modulation'   - Modulation type: 'sine', 'cosine', 'none' (default: 'sine')
%   'name'         - Configuration name (default: auto-generated)
%   'apply'        - Apply configuration immediately (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'carrier_freq', 100e3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'chip_rate', 100e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 400e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'modulation', 'sine', @(x) ismember(x, {'none', 'sine', 'cosine'}));
addParameter(p, 'name', '', @ischar);
addParameter(p, 'apply', true, @islogical);
addParameter(p, 'interactive', false, @islogical);

parse(p, varargin{:});
cfg = p.Results;

%% Interactive mode
if cfg.interactive || (nargin == 0)
    fprintf('=== SSTDR Custom Configuration ===\n');
    
    % Get modulation type
    fprintf('Available modulation types:\n');
    fprintf('  1. Sine wave modulation\n');
    fprintf('  2. Cosine wave modulation\n');
    fprintf('  3. Unmodulated (baseband)\n');
    mod_choice = input('Select modulation (1-3): ');
    
    switch mod_choice
        case 1
            cfg.modulation = 'sine';
        case 2
            cfg.modulation = 'cosine';
        case 3
            cfg.modulation = 'none';
        otherwise
            cfg.modulation = 'sine';
            fprintf('Invalid choice, using sine modulation\n');
    end
    
    % Get chip rate
    chip_rate_input = input(sprintf('Enter chip rate (Hz) [default: %.0f]: ', cfg.chip_rate));
    if ~isempty(chip_rate_input)
        cfg.chip_rate = chip_rate_input;
    end
    
    % Get carrier frequency
    if ~strcmp(cfg.modulation, 'none')
        carrier_input = input(sprintf('Enter carrier frequency (Hz) [default: %.0f]: ', cfg.carrier_freq));
        if ~isempty(carrier_input)
            cfg.carrier_freq = carrier_input;
        end
    else
        cfg.carrier_freq = 0;
    end
    
    % Get sampling frequency
    fs_input = input(sprintf('Enter sampling frequency (Hz) [default: %.0f]: ', cfg.fs));
    if ~isempty(fs_input)
        cfg.fs = fs_input;
    end
    
    name_input = input('Enter configuration name [auto-generate]: ', 's');
    if ~isempty(name_input)
        cfg.name = name_input;
    end
end

%% Validate frequency relationships
if ~strcmp(cfg.modulation, 'none')
    nyquist_freq = cfg.fs / 2;
    if cfg.carrier_freq > nyquist_freq
        warning('Carrier frequency (%.1f kHz) exceeds Nyquist frequency (%.1f kHz)', ...
                cfg.carrier_freq/1000, nyquist_freq/1000);
        fprintf('Consider increasing sampling frequency or reducing carrier frequency\n');
    end
    
    % Rule of thumb: sampling should be at least 4x chip rate
    if cfg.fs < 4 * cfg.chip_rate
        warning('Low sampling ratio (%.1fx). Recommend fs >= %.1f kHz for good quality', ...
                cfg.fs/cfg.chip_rate, 4*cfg.chip_rate/1000);
    end
    
    % Check carrier vs chip rate relationship
    if cfg.carrier_freq < cfg.chip_rate
        warning('Carrier frequency (%.1f kHz) is lower than chip rate (%.1f kHz)', ...
                cfg.carrier_freq/1000, cfg.chip_rate/1000);
        fprintf('This is unusual but may be intentional for your application\n');
    end
else
    % For unmodulated, just check sampling vs chip rate
    if cfg.fs < 4 * cfg.chip_rate
        warning('Low sampling ratio (%.1fx). Recommend fs >= %.1f kHz for good quality', ...
                cfg.fs/cfg.chip_rate, 4*cfg.chip_rate/1000);
    end
end

%% Generate configuration name
if isempty(cfg.name)
    if strcmp(cfg.modulation, 'none')
        cfg.name = sprintf('Custom_Unmod_%.0fkHz_Chip_%.0fkHz_Sampling', ...
                          cfg.chip_rate/1000, cfg.fs/1000);
    else
        cfg.name = sprintf('Custom_%s_%.0fkHz_Carrier_%.0fkHz_Chip_%.0fkHz_Sampling', ...
                          cfg.modulation, cfg.carrier_freq/1000, cfg.chip_rate/1000, cfg.fs/1000);
    end
end

%% Calculate optimal parameters based on frequencies
% Solver selection based on signal characteristics
if strcmp(cfg.modulation, 'none')
    solver = 'ode23t';   % Good for discontinuous signals
else
    solver = 'ode45';    % Good for smooth signals  
end

% Max step based on sampling frequency for numerical accuracy
max_step = 1/(10*cfg.fs);  % 10 simulation steps per sample period

% Stop time: enough for multiple PN periods
pn_length = 1023;  % 10-bit PN sequence
stop_time = pn_length / cfg.chip_rate;

%% Create configuration structure
config = struct( ...
    'name', cfg.name, ...
    'pn_config', struct( ...
        'modulation', cfg.modulation, ...
        'carrier_freq', cfg.carrier_freq, ...
        'chip_rate', cfg.chip_rate, ...
        'fs', cfg.fs, ...
        'pn_bits', 10, ...
        'polynomial', [10 3 0], ...
        'magnitude', 1 ...
    ), ...
    'correlation_config', struct( ...
        'method', 'freq', ...
        'peak_threshold', 0.1, ...
        'normalize', true, ...
        'plot_results', true ...
    ), ...
    'simulation_config', struct( ...
        'stop_time', stop_time, ...
        'solver', solver, ...
        'max_step', max_step ...
    ) ...
);

%% Display configuration summary
fprintf('\n=== Custom SSTDR Configuration ===\n');
fprintf('Name: %s\n', config.name);
fprintf('Modulation: %s\n', cfg.modulation);
fprintf('Frequencies:\n');
fprintf('  - Chip rate: %.1f kHz (%.3f μs per chip)\n', cfg.chip_rate/1000, 1e6/cfg.chip_rate);
if ~strcmp(cfg.modulation, 'none')
    fprintf('  - Carrier frequency: %.1f kHz\n', cfg.carrier_freq/1000);
    fprintf('  - Frequency ratio: %.1fx (carrier/chip)\n', cfg.carrier_freq/cfg.chip_rate);
end
fprintf('  - Sampling frequency: %.1f kHz\n', cfg.fs/1000);
fprintf('  - Sampling ratio: %.1fx (fs/chip_rate)\n', cfg.fs/cfg.chip_rate);
fprintf('Duration: %.3f ms\n', stop_time*1000);
fprintf('Solver: %s (max step: %.1f μs)\n', solver, max_step*1e6);

%% Apply configuration if requested
if cfg.apply
    fprintf('\nApplying configuration...\n');
    
    % Generate PN code with custom parameters
    pn_params = namedargs2cell(config.pn_config);
    pn_result = gen_pn_code(pn_params{:});
    
    % Store results
    config.pn_result = pn_result;
    
    % Export simulation parameters
    assignin('base', 'sim_stop_time', config.simulation_config.stop_time);
    assignin('base', 'sim_solver', config.simulation_config.solver);
    assignin('base', 'sim_max_step', config.simulation_config.max_step);
    assignin('base', 'sim_fs', config.pn_config.fs);
    assignin('base', 'sim_ts', 1/config.pn_config.fs);
    assignin('base', 'sim_chip_rate', config.pn_config.chip_rate);
    assignin('base', 'sim_decimation', pn_result.settings.interpFactor);
    assignin('base', 'sstdr_config', config);
    
    fprintf('✓ Configuration applied and variables exported\n');
    fprintf('✓ Ready for simulation!\n');
end

end 