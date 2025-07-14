function [config] = gen_pn_code(varargin)
%GEN_PN_CODE Generate PN sequence with configurable modulation for SSTDR
%
% Usage:
%   gen_pn_code()                           % Use default settings
%   gen_pn_code('modulation', 'sine')       % Sine wave modulation
%   gen_pn_code('modulation', 'none')       % Unmodulated
%   config = gen_pn_code(...)               % Return configuration
%
% Parameters:
%   'modulation'    - 'none', 'sine', 'cosine' (default: 'none')
%   'carrier_freq'  - Carrier frequency in Hz (default: 100e3)
%   'chip_rate'     - Chip rate in Hz (default: 125e3)
%   'fs'           - Sampling frequency in Hz (default: 1e6)
%   'interpFactor' - Interpolation factor (auto-calculated from fs/chip_rate if not specified)
%   'pn_bits'      - PN sequence bits (default: 10)
%   'polynomial'   - Primitive polynomial (default: [10 3 0])
%   'magnitude'    - Signal magnitude (default: 1)
%   'export_to_base' - Export variables to base workspace (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'modulation', 'none', @(x) ismember(x, {'none', 'sine', 'cosine'}));
addParameter(p, 'carrier_freq', 100e3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'chip_rate', 125e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 1e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'interpFactor', [], @(x) isnumeric(x) && x > 0);
addParameter(p, 'pn_bits', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'polynomial', [10 3 0], @(x) isnumeric(x));
addParameter(p, 'magnitude', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'export_to_base', true, @islogical);

parse(p, varargin{:});
cfg = p.Results;

%% Calculate interpolation factor from chip rate and sampling frequency
if isempty(cfg.interpFactor)
    % Calculate interpFactor to achieve desired chip rate
    cfg.interpFactor = round(cfg.fs / cfg.chip_rate);
    
    % Verify the calculation makes sense
    if cfg.interpFactor < 1
        cfg.interpFactor = 1;
        warning('Chip rate (%.1f kHz) higher than sampling rate (%.1f kHz). Setting interpFactor=1.', ...
                cfg.chip_rate/1000, cfg.fs/1000);
    end
    
    % Calculate actual chip rate achieved
    actual_chip_rate = cfg.fs / cfg.interpFactor;
    if abs(actual_chip_rate - cfg.chip_rate) / cfg.chip_rate > 0.01  % >1% error
        fprintf('Note: Requested chip rate %.1f kHz, achieved %.1f kHz (interpFactor=%d)\n', ...
                cfg.chip_rate/1000, actual_chip_rate/1000, cfg.interpFactor);
    end
    
    % Update config with actual chip rate
    cfg.chip_rate = actual_chip_rate;
else
    % If interpFactor was specified, calculate chip rate from it
    cfg.chip_rate = cfg.fs / cfg.interpFactor;
    fprintf('Using specified interpFactor=%d, chip rate = %.1f kHz\n', ...
            cfg.interpFactor, cfg.chip_rate/1000);
end

%% Generate PN sequence
fprintf('Generating %d-bit PN sequence...\n', cfg.pn_bits);

% Create PN sequence generator
% For polynomial [n a b c 0], the degree is n, so initial conditions need n elements
polynomial_degree = cfg.polynomial(1);
initial_conditions = [zeros(1, polynomial_degree-1) 1];

seq = comm.PNSequence( ...
    'Polynomial', cfg.polynomial, ...
    'SamplesPerFrame', 2^cfg.pn_bits - 1, ...
    'InitialConditions', initial_conditions);

% Generate one period of PN sequence
chips = 2*cfg.magnitude*seq() - cfg.magnitude;  % ±magnitude PN chips
pn_length = length(chips);

fprintf('PN sequence length: %d chips\n', pn_length);
fprintf('Chip rate: %.1f kHz (%.3f μs per chip)\n', cfg.chip_rate/1000, 1e6/cfg.chip_rate);

%% Create time vectors
Ts = 1/cfg.fs;
t_chip = (0:pn_length-1)' * Ts;  % Time vector for chip-rate sequence

% Interpolated time vector
code_interp = repelem(chips, cfg.interpFactor);
t_interp = (0:length(code_interp)-1)' * Ts;

%% Apply modulation
switch lower(cfg.modulation)
    case 'none'
        fprintf('Using unmodulated PN sequence\n');
        modulated_chips = chips;
        modulated_interp = code_interp;
        
    case 'sine'
        fprintf('Applying sine wave modulation at %.1f kHz\n', cfg.carrier_freq/1000);
        carrier_chip = sin(2*pi*cfg.carrier_freq*t_chip);
        carrier_interp = sin(2*pi*cfg.carrier_freq*t_interp);
        
        modulated_chips = chips .* carrier_chip;
        modulated_interp = code_interp .* carrier_interp;
        
    case 'cosine'
        fprintf('Applying cosine wave modulation at %.1f kHz\n', cfg.carrier_freq/1000);
        carrier_chip = cos(2*pi*cfg.carrier_freq*t_chip);
        carrier_interp = cos(2*pi*cfg.carrier_freq*t_interp);
        
        modulated_chips = chips .* carrier_chip;
        modulated_interp = code_interp .* carrier_interp;
end

%% Create TimeSeries objects for Simscape
codeTS_chip = timeseries(modulated_chips, t_chip);
codeTS_interp = timeseries(modulated_interp, t_interp);

%% Export to base workspace if requested
if cfg.export_to_base
    fprintf('Exporting variables to base workspace...\n');
    
    % Original variables for backward compatibility
    assignin('base', 'code_seq', chips);
    assignin('base', 'code_interp', code_interp);
    assignin('base', 'code_output', codeTS_chip);
    
    % New enhanced variables
    assignin('base', 'pn_chips', chips);
    assignin('base', 'pn_modulated', modulated_chips);
    assignin('base', 'pn_interp_modulated', modulated_interp);
    assignin('base', 'pn_timeseries_chip', codeTS_chip);
    assignin('base', 'pn_timeseries_interp', codeTS_interp);
    assignin('base', 'pn_config', cfg);
    
    fprintf('Variables exported:\n');
    fprintf('  - code_seq, code_interp, code_output (legacy)\n');
    fprintf('  - pn_chips, pn_modulated, pn_interp_modulated\n');
    fprintf('  - pn_timeseries_chip, pn_timeseries_interp\n');
    fprintf('  - pn_config (configuration struct)\n');
end

%% Create configuration output
config = struct();
config.settings = cfg;
config.pn_chips = chips;
config.pn_modulated = modulated_chips;
config.pn_interp_modulated = modulated_interp;
config.time_chip = t_chip;
config.time_interp = t_interp;
config.timeseries_chip = codeTS_chip;
config.timeseries_interp = codeTS_interp;
config.pn_length = pn_length;
config.total_duration = max(t_interp);
config.chip_rate = cfg.chip_rate;

fprintf('PN code generation complete!\n');
fprintf('Total duration: %.3f ms\n', config.total_duration*1000);
fprintf('Modulation: %s\n', cfg.modulation);
fprintf('Frequency summary:\n');
fprintf('  - Chip rate: %.1f kHz\n', cfg.chip_rate/1000);
if ~strcmp(cfg.modulation, 'none')
    fprintf('  - Carrier frequency: %.1f kHz\n', cfg.carrier_freq/1000);
end
fprintf('  - Sampling frequency: %.1f kHz\n', cfg.fs/1000);

end
