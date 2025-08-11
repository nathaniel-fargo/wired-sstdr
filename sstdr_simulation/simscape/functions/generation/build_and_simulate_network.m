function [results, model_info, sim_data] = build_and_simulate_network(network_config, simulation_config, varargin)
%BUILD_AND_SIMULATE_NETWORK Build network model and run SSTDR simulation
%
% This function combines the network configuration, model building, and
% simulation execution into a single workflow. It takes network and
% simulation configurations, builds a Simscape model, and runs the SSTDR
% simulation.
%
% Usage:
%   [results, model_info, sim_data] = build_and_simulate_network(network_config, simulation_config)
%   [results, model_info, sim_data] = build_and_simulate_network(network_config, simulation_config, 'save_model', true)
%
% Input:
%   network_config     - Network configuration structure from create_network_config
%   simulation_config  - Simulation configuration structure from create_simulation_config
%   varargin          - Optional name-value pairs:
%                       'model_name' - Name for the model (default: from network config)
%                       'save_model' - Save model to disk (default: false)
%                       'close_model' - Close model after simulation (default: true)
%                       'delete_model' - Delete model file after simulation (default: true)
%                       'connect_blocks' - Auto-connect blocks (default: true)
%                       'verbose' - Display progress messages (default: true)
%                       'save_data' - Save comprehensive data to .mat file (default: true)
%                       'data_folder' - Custom folder for data (default: auto-generated)
%
% Output:
%   results     - Simulation results structure with correlation analysis
%   model_info  - Model information structure with block details
%   sim_data    - Raw simulation data from Simulink

% Parse optional arguments
p = inputParser;
addRequired(p, 'network_config', @isstruct);
addRequired(p, 'simulation_config', @isstruct);
addParameter(p, 'model_name', '', @ischar);
addParameter(p, 'save_model', false, @islogical);
addParameter(p, 'close_model', true, @islogical);
addParameter(p, 'delete_model', true, @islogical);
addParameter(p, 'connect_blocks', true, @islogical);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_data', true, @islogical);
addParameter(p, 'data_folder', '', @ischar);
parse(p, network_config, simulation_config, varargin{:});
opts = p.Results;

% Use network name or generate model name
if isempty(opts.model_name)
    % Create unique model name based on network and timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    opts.model_name = sprintf('%s_%s', network_config.name, timestamp);
end

if opts.verbose
    fprintf('\n=== Build and Simulate Network Workflow ===\n');
    fprintf('Network: %s (%d segments, %.1f m)\n', ...
        network_config.name, network_config.num_segments, network_config.physical.total_length);
    fprintf('Simulation: %s\n', simulation_config.name);
    fprintf('Model: %s\n', opts.model_name);
end

% Initialize results structure
results = struct();
results.success = false;
results.timestamp = datestr(now);
results.network_config = network_config;
results.simulation_config = simulation_config;

try
    %% Step 1: Build the network model
    if opts.verbose
        fprintf('\n--- Step 1: Building Network Model ---\n');
    end
    
    % Save current directory and change to simscape root for model building
    original_dir = pwd;
    current_dir = fileparts(mfilename('fullpath'));
    simscape_root = fileparts(fileparts(current_dir));  % Go up from functions/generation
    cd(simscape_root);
    
    [model_name, model_info] = build_network_model(network_config, ...
        'model_name', opts.model_name, ...
        'save_model', opts.save_model, ...
        'close_after', false, ...  % Keep open for simulation
        'connect_blocks', opts.connect_blocks);
    
    % Restore original directory
    cd(original_dir);
    
    if opts.verbose
        fprintf('✓ Network model built successfully: %s\n', model_name);
    end
    
    %% Step 2: Run the simulation
    if opts.verbose
        fprintf('\n--- Step 2: Running SSTDR Simulation ---\n');
    end
    
    % Change to simscape root directory for simulation (same as model building)
    cd(simscape_root);
    
    [sim_results, sim_data] = run_simulation(simulation_config, model_name, ...
        'verbose', opts.verbose);
    
    % Restore original directory
    cd(original_dir);
    
    if opts.verbose
        fprintf('✓ Simulation completed successfully\n');
    end
    
    %% Step 3: Package results
    if opts.verbose
        fprintf('\n--- Step 3: Packaging Results ---\n');
    end
    
    % Combine all results
    results.success = true;
    results.model_name = model_name;
    results.simulation_results = sim_results;
    results.model_info = model_info;
    
    % Extract key simulation metrics
    if isfield(sim_results, 'sim_time')
        results.sim_time = sim_results.sim_time;
    end
    
    if isfield(sim_results, 'correlation_results')
        results.correlation_results = sim_results.correlation_results;
        
        % Extract peak information for easy access
        if isfield(sim_results.correlation_results, 'peak_times_freq')
            results.peak_times = sim_results.correlation_results.peak_times_freq;
            results.peak_values = sim_results.correlation_results.peaks_freq;
        elseif isfield(sim_results.correlation_results, 'peak_times_time')
            results.peak_times = sim_results.correlation_results.peak_times_time;
            results.peak_values = sim_results.correlation_results.peaks_time;
        end
    end
    
    % Network analysis summary
    results.network_analysis = struct();
    results.network_analysis.total_length = network_config.physical.total_length;
    results.network_analysis.num_faults = network_config.analysis.num_faults;
    results.network_analysis.fault_positions = network_config.analysis.fault_positions;
    results.network_analysis.max_fault_magnitude = network_config.analysis.max_fault_magnitude;
    
    if opts.verbose
        fprintf('✓ Results packaged successfully\n');
        
        % Display summary
        fprintf('\n=== Simulation Summary ===\n');
        fprintf('Network: %d segments, %d faults\n', ...
            network_config.num_segments, network_config.analysis.num_faults);
        fprintf('Simulation time: %.2f seconds\n', results.sim_time);
        
        if isfield(results, 'peak_times') && ~isempty(results.peak_times)
            fprintf('Detected peaks at: ');
            fprintf('%.3f ', results.peak_times * 1e6);
            fprintf('μs\n');
        end
        
        fprintf('Model: %s\n', model_name);
        fprintf('Status: ✓ Complete\n');
    end
    
    %% Step 3.5: Save comprehensive dataset
    if opts.save_data
        if opts.verbose
            fprintf('\n--- Step 3.5: Saving Dataset ---\n');
        end
        
        try
            [dataset_path, dataset_info] = save_sstdr_dataset(results, network_config, simulation_config, ...
                opts.data_folder, opts.verbose);
            
            results.dataset_path = dataset_path;
            results.dataset_info = dataset_info;
            
            if opts.verbose
                fprintf('✓ Dataset saved successfully\n');
                fprintf('  Path: %s\n', dataset_path);
            end
            
        catch ME
            if opts.verbose
                fprintf('⚠ Dataset saving failed: %s\n', ME.message);
            end
            % Don't fail the entire workflow if data saving fails
        end
    end
    
catch ME
    % Handle errors gracefully
    results.success = false;
    results.error = ME;
    
    % Restore original directory if it was changed
    if exist('original_dir', 'var')
        cd(original_dir);
    end
    
    if opts.verbose
        fprintf('✗ Workflow failed: %s\n', ME.message);
        fprintf('  Error in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    
    % Clean up model if it was created
    if exist('model_name', 'var') && bdIsLoaded(model_name)
        if opts.verbose
            fprintf('Cleaning up model: %s\n', model_name);
        end
        close_system(model_name, 0);  % Close without saving
        
        % Delete model file if requested
        if opts.delete_model
            delete_model_file(model_name);
            if opts.verbose
                fprintf('✓ Model file deleted during cleanup\n');
            end
        end
    end
    
    % Set empty outputs for consistency
    model_info = struct();
    sim_data = [];
    
    % Re-throw error if not in verbose mode (for debugging)
    if ~opts.verbose
        rethrow(ME);
    end
end

%% Step 4: Clean up
if opts.close_model && exist('model_name', 'var') && bdIsLoaded(model_name)
    if opts.verbose
        fprintf('\n--- Step 4: Cleanup ---\n');
        fprintf('Closing model: %s\n', model_name);
    end
    
    if opts.save_model
        % Ensure models directory exists
        models_dir = 'models';
        if ~exist(models_dir, 'dir')
            mkdir(models_dir);
        end
        
        % Save model to models subdirectory
        model_path = fullfile(models_dir, model_name);
        save_system(model_name, model_path);
        if opts.verbose
            fprintf('✓ Model saved: %s.slx\n', model_path);
        end
    end
    
    close_system(model_name, 0);  % Force close without saving changes
    if opts.verbose
        fprintf('✓ Model closed\n');
    end
    
    % Delete model file if requested
    if opts.delete_model
        delete_model_file(model_name);
        if opts.verbose
            fprintf('✓ Model file deleted\n');
        end
    end
end

if opts.verbose
    fprintf('\n=== Workflow Complete ===\n');
end

end

function delete_model_file(model_name)
%DELETE_MODEL_FILE Delete a Simulink model file from disk
%
% This function deletes both the .slx file and any associated .slxc files
% from both the root directory and models/ subdirectory.

    % Check both root directory and models/ subdirectory
    possible_paths = {
        [model_name '.slx'],
        [model_name '.slxc'],
        fullfile('models', [model_name '.slx']),
        fullfile('models', [model_name '.slxc'])
    };
    
    files_deleted = 0;
    
    for i = 1:length(possible_paths)
        if exist(possible_paths{i}, 'file')
            try
                delete(possible_paths{i});
                fprintf('  ✓ Deleted model file: %s\n', possible_paths{i});
                files_deleted = files_deleted + 1;
            catch ME
                fprintf('  ⚠ Could not delete %s: %s\n', possible_paths{i}, ME.message);
            end
        end
    end
    
    if files_deleted == 0
        fprintf('  ⚠ No model files found to delete for: %s\n', model_name);
    end
end 