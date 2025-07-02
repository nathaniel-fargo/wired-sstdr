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
%   'fs'           - Sampling frequency in Hz (default: 1e6)
%   'interpFactor' - Interpolation factor (default: 8)
%   'pn_bits'      - PN sequence bits (default: 10)
%   'polynomial'   - Primitive polynomial (default: [10 3 0])
%   'magnitude'    - Signal magnitude (default: 1)
%   'export_to_base' - Export variables to base workspace (default: true)

%% Parse input arguments
p = inputParser;
addParameter(p, 'modulation', 'none', @(x) ismember(x, {'none', 'sine', 'cosine'}));
addParameter(p, 'carrier_freq', 100e3, @(x) isnumeric(x) && x > 0);
addParameter(p, 'fs', 1e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'interpFactor', 8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'pn_bits', 10, @(x) isnumeric(x) && x > 0);
addParameter(p, 'polynomial', [10 3 0], @(x) isnumeric(x));
addParameter(p, 'magnitude', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'export_to_base', true, @islogical);

parse(p, varargin{:});
cfg = p.Results;

%% Generate PN sequence
fprintf('Generating %d-bit PN sequence...\n', cfg.pn_bits);

% Create PN sequence generator
seq = comm.PNSequence( ...
    'Polynomial', cfg.polynomial, ...
    'SamplesPerFrame', 2^cfg.pn_bits - 1, ...
    'InitialConditions', [zeros(1, cfg.pn_bits-1) 1]);

% Generate one period of PN sequence
chips = 2*cfg.magnitude*seq() - cfg.magnitude;  % ±magnitude PN chips
pn_length = length(chips);

fprintf('PN sequence length: %d chips\n', pn_length);

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

fprintf('PN code generation complete!\n');
fprintf('Total duration: %.3f ms\n', config.total_duration*1000);
fprintf('Modulation: %s\n', cfg.modulation);

end
