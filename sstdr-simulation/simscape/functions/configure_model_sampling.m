function configure_model_sampling(model_name, varargin)
%CONFIGURE_MODEL_SAMPLING Configure Simscape model sampling parameters
%
% Usage:
%   configure_model_sampling('my_model')           % Use current config
%   configure_model_sampling('my_model', 'default') % Load config first
%   configure_model_sampling('my_model', 'sine_fast', 'verbose', true)
%
% This function configures your Simscape model's sampling parameters
% to match the SSTDR configuration automatically.

%% Parse inputs
p = inputParser;
addRequired(p, 'model_name', @ischar);
addOptional(p, 'config_name', '', @ischar);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'auto_fix', true, @islogical);

parse(p, model_name, varargin{:});
opts = p.Results;

%% Load configuration if specified
if ~isempty(opts.config_name)
    if opts.verbose
        fprintf('Loading SSTDR configuration: %s\n', opts.config_name);
    end
    sstdr_config(opts.config_name);
end

%% Check if sampling parameters exist
if ~evalin('base', 'exist(''sim_ts'', ''var'')')
    error('Sampling parameters not found. Run sstdr_config() first.');
end

% Get sampling parameters
sim_fs = evalin('base', 'sim_fs');
sim_ts = evalin('base', 'sim_ts');
sim_decimation = evalin('base', 'sim_decimation');

if opts.verbose
    fprintf('Configuring model: %s\n', model_name);
    fprintf('  Sampling frequency: %.0f Hz\n', sim_fs);
    fprintf('  Sample time: %.2e s\n', sim_ts);
    fprintf('  Decimation factor: %d\n', sim_decimation);
end

%% Load model if not already loaded
if ~bdIsLoaded(model_name)
    if opts.verbose
        fprintf('Loading model: %s\n', model_name);
    end
    load_system(model_name);
end

%% Configure measurement blocks
measurement_blocks = find_system(model_name, 'BlockType', 'ToWorkspace');

for i = 1:length(measurement_blocks)
    block_name = measurement_blocks{i};
    
    % Check current sample time
    current_ts = get_param(block_name, 'SampleTime');
    
    if opts.verbose
        fprintf('Configuring measurement block: %s\n', block_name);
        fprintf('  Current sample time: %s\n', current_ts);
    end
    
    % Set sample time if needed
    if opts.auto_fix || strcmp(current_ts, '-1') || strcmp(current_ts, '0')
        set_param(block_name, 'SampleTime', 'sim_ts');
        if opts.verbose
            fprintf('  ✓ Set sample time to: sim_ts (%.2e s)\n', sim_ts);
        end
    else
        % Check if current setting is compatible
        try
            current_value = str2double(current_ts);
            if abs(current_value - sim_ts) > 1e-9
                if opts.auto_fix
                    set_param(block_name, 'SampleTime', 'sim_ts');
                    if opts.verbose
                        fprintf('  ✓ Fixed sample time mismatch\n');
                    end
                else
                    warning('Block %s has sample time %.2e, expected %.2e', ...
                            block_name, current_value, sim_ts);
                end
            else
                if opts.verbose
                    fprintf('  ✓ Sample time already correct\n');
                end
            end
        catch
            if opts.verbose
                fprintf('  → Sample time is variable: %s\n', current_ts);
            end
        end
    end
    
    % Configure other settings
    set_param(block_name, 'MaxDataPoints', 'inf');
    set_param(block_name, 'SaveFormat', 'Structure With Time');
end

%% Configure signal source blocks
source_blocks = find_system(model_name, 'BlockType', 'FromWorkspace');

for i = 1:length(source_blocks)
    block_name = source_blocks{i};
    
    if opts.verbose
        fprintf('Configuring source block: %s\n', block_name);
    end
    
    % Check if it's using SSTDR signals
    var_name = get_param(block_name, 'VariableName');
    if contains(var_name, 'pn_') || contains(var_name, 'code_')
        % Set interpolation method
        set_param(block_name, 'Interpolate', 'on');
        set_param(block_name, 'ZeroCross', 'on');
        
        if opts.verbose
            fprintf('  ✓ Configured for SSTDR signal: %s\n', var_name);
        end
    end
end

%% Configure solver settings
current_solver = get_param(model_name, 'Solver');
current_stop_time = str2double(get_param(model_name, 'StopTime'));

if evalin('base', 'exist(''sim_stop_time'', ''var'')')
    sim_stop_time = evalin('base', 'sim_stop_time');
    sim_solver = evalin('base', 'sim_solver');
    sim_max_step = evalin('base', 'sim_max_step');
    
    if opts.verbose
        fprintf('Solver configuration:\n');
        fprintf('  Current: %s, Stop: %.3f ms\n', current_solver, current_stop_time*1000);
        fprintf('  Recommended: %s, Stop: %.3f ms\n', sim_solver, sim_stop_time*1000);
    end
    
    if opts.auto_fix
        set_param(model_name, 'Solver', sim_solver);
        set_param(model_name, 'StopTime', num2str(sim_stop_time));
        set_param(model_name, 'MaxStep', num2str(sim_max_step));
        
        if opts.verbose
            fprintf('  ✓ Solver settings updated\n');
        end
    end
end

%% Summary
if opts.verbose
    fprintf('\n=== Model Configuration Summary ===\n');
    fprintf('Model: %s\n', model_name);
    fprintf('Measurement blocks configured: %d\n', length(measurement_blocks));
    fprintf('Source blocks configured: %d\n', length(source_blocks));
    fprintf('Sampling rate: %.0f Hz (%.2e s)\n', sim_fs, sim_ts);
    
    if opts.auto_fix
        fprintf('Status: ✓ All settings automatically configured\n');
    else
        fprintf('Status: → Manual configuration required\n');
    end
    
    fprintf('Ready for simulation!\n');
end

end 