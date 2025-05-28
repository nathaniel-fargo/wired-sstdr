% display_network_results.m
%
% Combines network diagram and LiveWire results into a single display.
% Can also save the combined plot to a file if an output folder is provided.
%
% Author: Nathaniel Fargo
% Date: 2024-07-30
% Org: U of U WIRED
%
% Usage:
%   display_network_results(networks_csv_file, livewire_results_folder);
%   display_network_results(networks_csv_file, livewire_results_folder, {'T0', 'T1'});
%   display_network_results(networks_csv_file, livewire_results_folder, {}, 'CombinedPlots/');
%   display_network_results(networks_csv_file, livewire_results_folder, {'T0'}, 'CombinedPlots/');

function display_network_results(networks_csv_file, livewire_results_folder, specific_network_ids, output_save_folder)

    close all;

    % --- Argument Handling ---
    if nargin < 2
        error('Usage: display_network_results(networks_csv_file, livewire_results_folder, [specific_network_ids], [output_save_folder])');
    end
    if nargin < 3
        specific_network_ids = []; % Process all networks by default
    end
    if nargin < 4
        output_save_folder = ''; % Display plots, don''t save by default
    end

    % --- Validate Inputs ---
    if ~isfile(networks_csv_file)
        error('Networks CSV file not found: %s', networks_csv_file);
    end
    if ~isfolder(livewire_results_folder)
        error('LiveWire results folder not found: %s', livewire_results_folder);
    end
    if ~isempty(output_save_folder) && ~isfolder(output_save_folder)
        fprintf('Output save folder %s does not exist. Creating it...\n', output_save_folder);
        mkdir(output_save_folder);
    end

    % --- Read Networks Data ---
    try
        netTable = readtable(networks_csv_file, 'Delimiter', ',', 'TextType', 'string');
    catch ME
        error('Failed to read networks CSV file: %s\n%s', networks_csv_file, ME.message);
    end
    
    networkIDsFromFile = netTable{:,1}; % Assuming ID is the first column

    % --- Filter Networks ---
    if ~isempty(specific_network_ids)
        if ischar(specific_network_ids) % Convert single char ID to cell
            specific_network_ids = {specific_network_ids};
        end
        if ~iscellstr(specific_network_ids) && ~isstring(specific_network_ids) %iscellstr for older MATLAB
             error('specific_network_ids must be a cell array of strings or a string array.');
        end
        specific_network_ids = string(specific_network_ids); % Ensure string array for consistent comparison
        
        [lia, ~] = ismember(lower(specific_network_ids), lower(networkIDsFromFile));
        if ~all(lia)
            warning('Some specified network IDs were not found in %s: %s', ...
                    networks_csv_file, strjoin(specific_network_ids(~lia), ', '));
        end
        networkIDsToProcess = specific_network_ids(lia);
        if isempty(networkIDsToProcess)
            disp('No valid network IDs to process after filtering.');
            return;
        end
    else
        networkIDsToProcess = networkIDsFromFile;
    end

    if isempty(networkIDsToProcess)
        disp('No networks to process.');
        return;
    end

    disp(['Processing ' num2str(length(networkIDsToProcess)) ' network(s)...']);

    % --- Loop Through Networks and Plot ---
    for i = 1:length(networkIDsToProcess)
        current_network_id_str = char(networkIDsToProcess(i)); % Use char for path construction and titles
        
        fprintf('\nProcessing Network ID: %s\n', current_network_id_str);

        % Construct path to the LiveWire data CSV for the current network
        livewire_csv_filename = [current_network_id_str, '.csv'];
        livewire_csv_filepath = fullfile(livewire_results_folder, livewire_csv_filename);

        if ~isfile(livewire_csv_filepath)
            warning('LiveWire results CSV file not found for network %s: %s. Skipping this network.', ...
                    current_network_id_str, livewire_csv_filepath);
            continue;
        end

        % --- Create Figure and Subplots ---
        fig = figure('Name', ['Network Analysis: ' current_network_id_str], 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);
        
        ax_network_diagram = subplot(2, 1, 1); 
        ax_livewire_data = subplot(2, 1, 2);

        % --- Plot Network Diagram ---
        fprintf('  Drawing network diagram for %s...\n', current_network_id_str);
        try
            draw_network(networks_csv_file, {current_network_id_str}, '', ax_network_diagram);
            title(ax_network_diagram, ['Network Diagram: ' current_network_id_str]);
        catch ME_draw
            warning('Error drawing network %s: %s. Diagram might be incomplete.', current_network_id_str, ME_draw.message);
            title(ax_network_diagram, ['Network Diagram: ' current_network_id_str ' (Error)']);
            text(ax_network_diagram, 0.5, 0.5, 'Error drawing network.', 'HorizontalAlignment', 'center');
        end

        % --- Plot LiveWire Data ---
        fprintf('  Plotting LiveWire data for %s from %s...\n', current_network_id_str, livewire_csv_filepath);
        try
            plot_livewire_csv(livewire_csv_filepath, {}, '', ax_livewire_data);
        catch ME_plot
            warning('Error plotting LiveWire data for %s: %s. Plot might be incomplete.', current_network_id_str, ME_plot.message);
            title(ax_livewire_data, ['LiveWire Data: ' current_network_id_str ' (Error)']);
            text(ax_livewire_data, 0.5, 0.5, 'Error plotting LiveWire data.', 'HorizontalAlignment', 'center');
        end
        
        % Add an overall super title to the figure
        % try
        %     sgtitle(fig, ['Combined Analysis for Network: ' current_network_id_str]);
        % catch
        %     annotation(fig, 'textbox', [0, 0.95, 1, 0.05], ...
        %                'String', ['Combined Analysis for Network: ' current_network_id_str], ...
        %                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        % end

        % --- Save Figure (Optional) ---
        if ~isempty(output_save_folder)
            output_filename = fullfile(output_save_folder, ['Combined_' current_network_id_str '.png']);
            try
                saveas(fig, output_filename);
                fprintf('  Combined plot saved to: %s\n', output_filename);
                close(fig); % Close after saving
            catch ME_save
                warning('Failed to save combined plot %s: %s', output_filename, ME_save.message);
                if isempty(output_save_folder) && ishandle(fig) 
                   % if saving was not requested, it should remain open
                elseif ishandle(fig)
                   close(fig); 
                end
            end
        else
            disp('  Displaying combined plot. Close figure to continue...');
        end
    end

    disp('\nFinished processing all selected networks.');
end
