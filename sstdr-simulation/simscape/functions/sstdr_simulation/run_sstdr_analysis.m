function run_sstdr_analysis(config_name, sim_model, varargin)
%RUN_SSTDR_ANALYSIS Complete SSTDR simulation and analysis workflow
%
% Usage:
%   run_sstdr_analysis()                           % Default config, no simulation
%   run_sstdr_analysis('sine_fast')                % Use sine_fast config
%   run_sstdr_analysis('default', 'my_model')      % Run simulation with model
%   run_sstdr_analysis(config, model, 'method', 'both')  % Specify correlation method
%
% Parameters:
%   config_name     - Configuration name or 'skip' to use existing workspace
%   sim_model       - Simscape model name (optional, '' to skip simulation)  
%   'method'        - Correlation method override ('time', 'freq', 'both')
%   'auto_run'      - Automatically run simulation (default: true)
%   'save_results'  - Save results to file (default: false)

%% Setup paths for organized folder structure
current_dir = fileparts(mfilename('fullpath'));
simscape_root = fileparts(fileparts(current_dir));  % Go up two levels from functions/sstdr_simulation
addpath(fullfile(simscape_root, 'functions', 'network_generation'));
addpath(fullfile(simscape_root, 'functions', 'sstdr_simulation'));
addpath(fullfile(simscape_root, 'config'));

%% Parse inputs
if nargin < 1 || isempty(config_name)
    config_name = 'default';
end

if nargin < 2
    sim_model = 'sstdr_basic';
end

p = inputParser;
addParameter(p, 'method', '', @ischar);
addParameter(p, 'auto_run', true, @islogical);
addParameter(p, 'save_results', false, @islogical);
addParameter(p, 'plot_results', [], @(x) isempty(x) || islogical(x));

parse(p, varargin{:});
opts = p.Results;

fprintf('=== SSTDR Analysis Workflow ===\n');
fprintf('Started at: %s\n', datestr(now));

%% Step 1: Load Configuration
if ~strcmp(config_name, 'skip')
    fprintf('\n--- Step 1: Loading Configuration ---\n');
    try
        config = sstdr_config(config_name);
        fprintf('✓ Configuration loaded successfully\n');
    catch ME
        fprintf('✗ Configuration failed: %s\n', ME.message);
        return;
    end
else
    fprintf('\n--- Step 1: Using Existing Workspace ---\n');
    if evalin('base', 'exist(''sstdr_config'', ''var'')')
        config = evalin('base', 'sstdr_config');
        fprintf('✓ Using existing configuration: %s\n', config.name);
    else
        fprintf('✗ No existing configuration found\n');
        return;
    end
end

%% Step 2: Run Simulation (if model provided)
if ~isempty(sim_model) && opts.auto_run
    fprintf('\n--- Step 2: Running Simscape Simulation ---\n');
    try
        % Check if model exists
        if ~exist(sim_model, 'file')
            if exist([sim_model '.slx'], 'file')
                sim_model = [sim_model '.slx'];
            else
                error('Model file not found: %s', sim_model);
            end
        end
        
        fprintf('Loading model: %s\n', sim_model);
        load_system(sim_model);
        
        % Apply simulation parameters
        set_param(sim_model, 'StopTime', num2str(config.simulation_config.stop_time));
        set_param(sim_model, 'Solver', config.simulation_config.solver);
        set_param(sim_model, 'MaxStep', num2str(config.simulation_config.max_step));
        
        fprintf('Running simulation (%.1f ms duration)...\n', config.simulation_config.stop_time*1000);
        
        % Run simulation
        tic;
        sim_out = sim(sim_model);
        sim_time = toc;
        
        % Store results in base workspace
        assignin('base', 'out', sim_out);
        
        fprintf('✓ Simulation completed in %.2f seconds\n', sim_time);
        
    catch ME
        fprintf('✗ Simulation failed: %s\n', ME.message);
        fprintf('  Continuing with existing simulation data...\n');
    end
elseif ~isempty(sim_model)
    fprintf('\n--- Step 2: Manual Simulation Required ---\n');
    fprintf('Please manually run simulation with model: %s\n', sim_model);
    fprintf('Recommended settings:\n');
    fprintf('  Stop Time: %.6f\n', config.simulation_config.stop_time);
    fprintf('  Solver: %s\n', config.simulation_config.solver);
    fprintf('  Max Step: %.2e\n', config.simulation_config.max_step);
    input('Press Enter when simulation is complete...');
else
    fprintf('\n--- Step 2: Using Existing Simulation Data ---\n');
    if evalin('base', 'exist(''out'', ''var'')')
        fprintf('✓ Found existing simulation data\n');
    else
        fprintf('⚠ No simulation data found. Analysis may fail.\n');
    end
end

%% Step 3: Cross-Correlation Analysis
fprintf('\n--- Step 3: Cross-Correlation Analysis ---\n');

% Determine correlation method
if ~isempty(opts.method)
    correlation_method = opts.method;
else
    correlation_method = config.correlation_config.method;
end

try
    % Get correlation parameters from config
    corr_params = namedargs2cell(config.correlation_config);
    
    % Override method if specified
    if ~isempty(opts.method)
        corr_params{find(strcmp(corr_params, 'method'))+1} = opts.method;
    end
    
    % Override plot_results if specified
    if ~isempty(opts.plot_results)
        plot_idx = find(strcmp(corr_params, 'plot_results'));
        if ~isempty(plot_idx)
            corr_params{plot_idx+1} = opts.plot_results;
        else
            corr_params{end+1} = 'plot_results';
            corr_params{end+1} = opts.plot_results;
        end
    end
    
    fprintf('Running %s domain correlation...\n', correlation_method);
    
    % Run correlation analysis
    tic;
    results = cross_correlate(corr_params{:});
    analysis_time = toc;
    
    fprintf('✓ Correlation analysis completed in %.2f seconds\n', analysis_time);
    
    % Display peak information
    if isfield(results, 'peak_times_time') && ~isempty(results.peak_times_time)
        fprintf('Time domain peaks at: ');
        fprintf('%.3f ', results.peak_times_time * 1e6);
        fprintf('µs\n');
    end
    
    if isfield(results, 'peak_times_freq') && ~isempty(results.peak_times_freq)
        fprintf('Frequency domain peaks at: ');
        fprintf('%.3f ', results.peak_times_freq * 1e6);
        fprintf('µs\n');
    end
    
catch ME
    fprintf('✗ Correlation analysis failed: %s\n', ME.message);
    results = [];
end

%% Step 4: Save Results (if requested)
if opts.save_results && ~isempty(results)
    fprintf('\n--- Step 4: Saving Results ---\n');
    try
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('sstdr_results_%s_%s.mat', config_name, timestamp);
        
        save(filename, 'config', 'results', 'timestamp');
        fprintf('✓ Results saved to: %s\n', filename);
        
    catch ME
        fprintf('✗ Failed to save results: %s\n', ME.message);
    end
end

%% Summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Configuration: %s\n', config.name);
fprintf('PN Length: %d chips\n', config.pn_result.pn_length);
fprintf('Modulation: %s\n', config.pn_config.modulation);
fprintf('Correlation Method: %s\n', correlation_method);

if ~isempty(results)
    if isfield(results, 'speedup')
        fprintf('Performance: %.1fx speedup with frequency domain\n', results.speedup);
    end
    
    if isfield(results, 'method_correlation')
        fprintf('Method Agreement: %.4f correlation\n', results.method_correlation);
    end
    
    fprintf('Status: ✓ Analysis completed successfully\n');
else
    fprintf('Status: ✗ Analysis failed\n');
end

fprintf('Completed at: %s\n', datestr(now));

end 