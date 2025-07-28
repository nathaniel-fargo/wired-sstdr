% test_measurements_parsing.m
%
% Simple test to verify measurements.csv parsing is working correctly
%
% Author: Nathaniel Fargo  
% Date: 2025-01-23

clear; close all; clc;

fprintf('Testing measurements.csv parsing...\n\n');

try
    % Test reading measurements.csv with the fixed approach
    measurements = readtable('measurements.csv', 'Delimiter', ';', 'VariableNamingRule', 'preserve');
    fprintf('✓ Successfully read measurements.csv\n');
    
    % Display table info
    fprintf('Table dimensions: %d rows × %d columns\n', height(measurements), width(measurements));
    fprintf('Column names: %s\n', strjoin(measurements.Properties.VariableNames, ', '));
    
    % Convert to cell arrays for easier handling
    ids = table2cell(measurements(:, 1));
    heights = table2cell(measurements(:, 2));
    terms = table2cell(measurements(:, 3));
    
    fprintf('\n--- Measurement Data ---\n');
    for i = 1:length(ids)
        id = ids{i};
        height_val = heights{i};
        term_val = terms{i};
        
        % Convert height to string if it's numeric or other type
        if isnumeric(height_val)
            height_str = num2str(height_val);
        elseif ischar(height_val) || isstring(height_val)
            height_str = char(height_val);
        else
            height_str = 'NA';
        end
        
        % Convert termination to string
        if ischar(term_val) || isstring(term_val)
            term_str = char(term_val);
        else
            term_str = 'Unknown';
        end
        
        fprintf('%s: %s inches, Termination: %s\n', id, height_str, term_str);
    end
    
    fprintf('\n✓ Successfully parsed all data!\n');
    
    % Test a quick plot to make sure everything works
    fprintf('\nTesting plot generation...\n');
    addpath('scripts');
    
    if isfolder('06-16-2025/LiveWire/CSV/')
        fig = plot_livewire_folder_time('06-16-2025/LiveWire/CSV/', '48 MHz');
        title('Test Plot - Measurements.csv Parsing Fixed');
        saveas(fig, 'measurements_parsing_test.png');
        close(fig);
        fprintf('✓ Test plot saved as measurements_parsing_test.png\n');
    else
        fprintf('⚠ Data folder not found, skipping plot test\n');
    end
    
catch ME
    fprintf('✗ Error: %s\n', char(ME.message));
    fprintf('Full error details:\n');
    disp(ME);
end

fprintf('\nTest completed!\n'); 