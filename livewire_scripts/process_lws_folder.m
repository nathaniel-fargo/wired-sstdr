% process_lws_folder.m
%
% Processes all .lws files in a specified input directory.
% 1. Converts .lws files to .csv files and stores them in a 'CSV_Files' subdirectory.
% 2. Generates plots from these .csv files and stores them in a 'Plot_Files' subdirectory.
%
% Author: Nathaniel Fargo
% Date: 2025-05-21
% Org: U of U WIRED
%
% Usage:
%   process_lws_folder('path/to/your/lws_folder');

function process_lws_folder(inputDir)

    if nargin < 1
        error('Usage: process_lws_folder(inputDir)');
    end

    if ~exist(inputDir, 'dir')
        error('Input LWS directory not found: %s', inputDir);
    end

    % Define subdirectory names
    csvSubDirName = 'CSV';
    plotsSubDirName = 'Plots';

    % Rename the input directory to 'LWS' and update inputDir
    parentDir = fileparts(inputDir);
    newLwsDir = fullfile(parentDir, 'LWS');
    if ~exist(newLwsDir, 'dir')
        fprintf('Renaming input directory %s to %s...\n', inputDir, newLwsDir);
        movefile(inputDir, newLwsDir);
    else
        fprintf('Directory %s already exists. Cannot rename %s.\n', newLwsDir, inputDir);
        error('Renaming failed: target directory already exists.');
    end
    inputDir = newLwsDir;

    % Create full paths for subdirectories in parent directory
    outputCsvDir = fullfile(parentDir, csvSubDirName);
    outputPlotsDir = fullfile(parentDir, plotsSubDirName);

    % Create subdirectories if they don't exist
    if ~exist(outputCsvDir, 'dir')
        fprintf('CSV output directory %s does not exist. Creating it...\n', outputCsvDir);
        mkdir(outputCsvDir);
    else
        fprintf('CSV output directory %s already exists.\n', outputCsvDir);
    end

    if ~exist(outputPlotsDir, 'dir')
        fprintf('Plots output directory %s does not exist. Creating it...\n', outputPlotsDir);
        mkdir(outputPlotsDir);
    else
        fprintf('Plots output directory %s already exists.\n', outputPlotsDir);
    end

    % Step 1: Convert .lws files to .csv
    fprintf('\nStep 1: Converting .lws files to .csv format...\n');
    fprintf('Input LWS directory: %s\n', inputDir);
    fprintf('Output CSV directory: %s\n', outputCsvDir);
    
    try
        lws_to_csv(inputDir, outputCsvDir);
        fprintf('Successfully converted .lws files to .csv files.\n');
    catch ME
        fprintf('ERROR during LWS to CSV conversion: %s\n', ME.message);
        % disp(ME.getReport());
        fprintf('Aborting further processing.\n');
        return;
    end

    % Step 2: Generate plots from .csv files
    fprintf('\nStep 2: Generating plots from .csv files...\n');
    fprintf('Input CSV directory: %s\n', outputCsvDir);
    fprintf('Output Plots directory: %s\n', outputPlotsDir);

    try
        batch_plot_livewire_csv(outputCsvDir, outputPlotsDir);
        fprintf('Successfully generated plots.\n');
    catch ME
        fprintf('ERROR during plot generation: %s\n', ME.message);
        % disp(ME.getReport());
    end

    fprintf('\nProcessing complete for folder: %s\n', inputDir);

end 