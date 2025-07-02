function [config] = sstdr_custom_config(varargin)
%SSTDR_CUSTOM_CONFIG Create custom SSTDR configuration with specified frequencies
%
% Usage:
%   sstdr_custom_config()                                    % Interactive mode
%   sstdr_custom_config('carrier_freq', 500e3)               % 500 kHz carrier
%   sstdr_custom_config('fs', 5e6)                          % 5 MHz sampling
%   sstdr_custom_config('carrier_freq', 1e6, 'fs', 10e6)    % Both frequencies
%   config = sstdr_custom_config(...)                       % Return config
%
% Parameters:
%   'carrier_freq'  - Carrier frequency in Hz (default: 100e3)
%   'fs'           - Sampling frequency in Hz (default: 2e6)  
%   'modulation'   - Modulation type: 'sine', 'cosine', 'none' (default: 'sine')
%   'name'         - Configuration name (default: auto-generated)
%   'apply'        - Apply configuration immediately (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'carrier_freq', 100e3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'fs', 2e6, @(x) isnumeric(x) && x > 0);
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
    
    % Get frequencies
    if ~strcmp(cfg.modulation, 'none')
        carrier_input = input(sprintf('Enter carrier frequency (Hz) [default: %.0f]: ', cfg.carrier_freq));
        if ~isempty(carrier_input)
            cfg.carrier_freq = carrier_input;
        end
    else
        cfg.carrier_freq = 0;
    end
    
    fs_input = input(sprintf('Enter sampling frequency (Hz) [default: %.0f]: ', cfg.fs));
    if ~isempty(fs_input)
        cfg.fs = fs_input;
    end
    
    name_input = input('Enter configuration name [auto-generate]: ', 's');
    if ~isempty(name_input)
        cfg.name = name_input;
    end
end

%% Validate frequency relationship
if ~strcmp(cfg.modulation, 'none')
    nyquist_freq = cfg.fs / 2;
    if cfg.carrier_freq > nyquist_freq
        warning('Carrier frequency (%.1f kHz) exceeds Nyquist frequency (%.1f kHz)', ...
                cfg.carrier_freq/1000, nyquist_freq/1000);
        fprintf('Consider increasing sampling frequency or reducing carrier frequency\n');
    end
    
    % Rule of thumb: sampling should be at least 10x carrier for good quality
    if cfg.fs < 10 * cfg.carrier_freq
        warning('Low sampling ratio (%.1fx). Consider fs >= %.1f MHz for better quality', ...
                cfg.fs/cfg.carrier_freq, 10*cfg.carrier_freq/1e6);
    end
end

%% Generate configuration name
if isempty(cfg.name)
    if strcmp(cfg.modulation, 'none')
        cfg.name = sprintf('Custom_Unmod_%.0fkHz_Sampling', cfg.fs/1000);
    else
        cfg.name = sprintf('Custom_%s_%.0fkHz_Carrier_%.0fkHz_Sampling', ...
                          cfg.modulation, cfg.carrier_freq/1000, cfg.fs/1000);
    end
end

%% Calculate optimal parameters
% Interpolation factor: balance between resolution and computation
if cfg.fs >= 5e6
    interp_factor = 16;  % High resolution
elseif cfg.fs >= 2e6
    interp_factor = 8;   % Standard resolution
else
    interp_factor = 4;   % Fast computation
end

% Solver selection based on signal characteristics
if strcmp(cfg.modulation, 'none')
    solver = 'ode23t';   % Good for discontinuous signals
    max_step = 1e-6;
else
    solver = 'ode45';    % Good for smooth signals  
    max_step = 1/(10*cfg.carrier_freq);  % 10 samples per carrier period
end

% Stop time: enough for multiple PN periods
pn_length = 1023;  % 10-bit PN sequence
stop_time = (pn_length * interp_factor) / cfg.fs;

%% Create configuration structure
config = struct( ...
    'name', cfg.name, ...
    'pn_config', struct( ...
        'modulation', cfg.modulation, ...
        'carrier_freq', cfg.carrier_freq, ...
        'fs', cfg.fs, ...
        'interpFactor', interp_factor, ...
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
if ~strcmp(cfg.modulation, 'none')
    fprintf('Carrier frequency: %.1f kHz\n', cfg.carrier_freq/1000);
    fprintf('Sampling ratio: %.1fx (fs/fc)\n', cfg.fs/cfg.carrier_freq);
end
fprintf('Sampling frequency: %.1f kHz\n', cfg.fs/1000);
fprintf('Interpolation factor: %d\n', interp_factor);
fprintf('Duration: %.3f ms\n', stop_time*1000);
fprintf('Solver: %s (max step: %.1f µs)\n', solver, max_step*1e6);

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
    assignin('base', 'sim_decimation', config.pn_config.interpFactor);
    assignin('base', 'sstdr_config', config);
    
    fprintf('✓ Configuration applied and variables exported\n');
    fprintf('✓ Ready for simulation!\n');
end

end 