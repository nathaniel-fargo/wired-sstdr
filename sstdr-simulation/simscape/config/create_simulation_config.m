function [config] = create_simulation_config(varargin)
%CREATE_SIMULATION_CONFIG Create custom SSTDR configuration with specified frequencies
%
% Usage:
%   create_simulation_config()                                    % Interactive mode
%   create_simulation_config('carrier_freq', 500e3)               % 500 kHz carrier
%   create_simulation_config('chip_rate', 200e3)                  % 200 kHz chip rate
%   create_simulation_config('fs', 5e6)                          % 5 MHz sampling
%   create_simulation_config('carrier_freq', 1e6, 'chip_rate', 1e6, 'fs', 4e6) % All frequencies
%   config = create_simulation_config(...)                       % Return config
%
% Parameters:
%   'carrier_freq'  - Carrier frequency in Hz (default: 100e3)
%   'chip_rate'     - Chip rate in Hz (default: 100e3)
%   'fs'           - Sampling frequency in Hz (default: 400e3)  
%   'modulation'   - Modulation type: 'sine', 'cosine', 'none' (default: 'sine')
%   'pn_bits'      - PN sequence bits (default: 10, gives 2^10-1 = 1023 chips)
%   'duration'     - Simulation duration in seconds (default: auto from PN length)
%   'name'         - Configuration name (default: auto-generated)
%   'apply'        - Apply configuration immediately (default: true)
%   'method'       - Cross-correlation method: 'time', 'freq', 'both' (default: 'freq')
%   'plot_results' - Show correlation plots (default: true)
%   'positive_only'- Return only positive time lags (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'carrier_freq', 100e3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'chip_rate', 100e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 400e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'modulation', 'sine', @(x) ismember(x, {'none', 'sine', 'cosine'}));
addParameter(p, 'pn_bits', 10, @(x) isnumeric(x) && x > 0 && x <= 20);
addParameter(p, 'duration', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'name', '', @ischar);
addParameter(p, 'apply', true, @islogical);
addParameter(p, 'interactive', false, @islogical);
addParameter(p, 'method', 'freq', @(x) ismember(x, {'time', 'freq', 'both'}));
addParameter(p, 'plot_results', true, @islogical);
addParameter(p, 'positive_only', true, @islogical);

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
    
    % Get PN sequence bits
    pn_bits_input = input(sprintf('Enter PN sequence bits [default: %d, gives %d chips]: ', cfg.pn_bits, 2^cfg.pn_bits - 1));
    if ~isempty(pn_bits_input)
        cfg.pn_bits = pn_bits_input;
    end
    
    % Get simulation duration
    if isempty(cfg.duration)
        default_duration = (2^cfg.pn_bits - 1) / cfg.chip_rate;
        duration_input = input(sprintf('Enter simulation duration (s) [default: %.2e s]: ', default_duration));
        if ~isempty(duration_input)
            cfg.duration = duration_input;
        else
            cfg.duration = default_duration;
        end
    end
    
    name_input = input('Enter configuration name [auto-generate]: ', 's');
    if ~isempty(name_input)
        cfg.name = name_input;
    end
    
    % Get correlation method
    fprintf('\nCorrelation Analysis Options:\n');
    fprintf('Available correlation methods:\n');
    fprintf('  1. Frequency domain (recommended)\n');
    fprintf('  2. Time domain\n');
    fprintf('  3. Both (comparison)\n');
    method_choice = input('Select correlation method (1-3): ');
    
    switch method_choice
        case 1
            cfg.method = 'freq';
        case 2
            cfg.method = 'time';
        case 3
            cfg.method = 'both';
        otherwise
            cfg.method = 'freq';
            fprintf('Invalid choice, using frequency domain\n');
    end
    
    % Get plot results option
    plot_input = input('Show correlation plots? (y/n) [y]: ', 's');
    if isempty(plot_input) || strcmpi(plot_input, 'y')
        cfg.plot_results = true;
    else
        cfg.plot_results = false;
    end
    
    % Get positive_only option
    positive_input = input('Return only positive time lags? (y/n) [y]: ', 's');
    if isempty(positive_input) || strcmpi(positive_input, 'y')
        cfg.positive_only = true;
    else
        cfg.positive_only = false;
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

% Stop time: use configured duration or calculate from PN length
pn_length = 2^cfg.pn_bits - 1;  % Configurable PN sequence length
if isempty(cfg.duration)
    stop_time = pn_length / cfg.chip_rate;
else
    stop_time = cfg.duration;
end

% Select appropriate polynomial based on pn_bits
polynomial = get_pn_polynomial(cfg.pn_bits);

%% Create configuration structure
config = struct( ...
    'name', cfg.name, ...
    'pn_config', struct( ...
        'modulation', cfg.modulation, ...
        'carrier_freq', cfg.carrier_freq, ...
        'chip_rate', cfg.chip_rate, ...
        'fs', cfg.fs, ...
        'pn_bits', cfg.pn_bits, ...
        'polynomial', polynomial, ...
        'magnitude', 1 ...
    ), ...
    'correlation_config', struct( ...
        'method', cfg.method, ...
        'peak_threshold', 0.1, ...
        'normalize', true, ...
        'plot_results', cfg.plot_results, ...
        'positive_only', cfg.positive_only ...
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
if cfg.plot_results
    plot_str = ', plots enabled';
else
    plot_str = ', plots disabled';
end
if cfg.positive_only
    lag_str = ', positive lags only';
else
    lag_str = ', full lags';
end
fprintf('Correlation: %s domain%s%s\n', cfg.method, plot_str, lag_str);

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

function polynomial = get_pn_polynomial(pn_bits)
%GET_PN_POLYNOMIAL Get appropriate primitive polynomial for given PN bits
%
% Common primitive polynomials for different PN sequence lengths

switch pn_bits
    case 3
        polynomial = [3 2 0];
    case 4
        polynomial = [4 3 0];
    case 5
        polynomial = [5 3 0];
    case 6
        polynomial = [6 5 0];
    case 7
        polynomial = [7 6 0];
    case 8
        polynomial = [8 6 5 4 0];
    case 9
        polynomial = [9 5 0];
    case 10
        polynomial = [10 3 0];
    case 11
        polynomial = [11 2 0];
    case 12
        polynomial = [12 6 4 1 0];
    case 13
        polynomial = [13 4 3 1 0];
    case 14
        polynomial = [14 5 3 1 0];
    case 15
        polynomial = [15 14 0];
    case 16
        polynomial = [16 5 3 2 0];
    case 17
        polynomial = [17 3 0];
    case 18
        polynomial = [18 7 0];
    case 19
        polynomial = [19 5 2 1 0];
    case 20
        polynomial = [20 3 0];
    otherwise
        % Default to 10-bit polynomial
        polynomial = [10 3 0];
        warning('Unsupported pn_bits value %d, using default 10-bit polynomial', pn_bits);
end
end 