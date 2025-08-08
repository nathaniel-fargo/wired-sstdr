% process_lws_folder.m
%
% Processes all .lws files in a specified input directory.
% 1. Converts .lws files to .csv files and stores them in a 'CSV' subdirectory.
% 2. Generates plots from these .csv files and stores them in a 'Plots' subdirectory.
%
% Author: Nathaniel Fargo
% Date: 2025-05-21
% Org: U of U WIRED
%
% Developed partially with AI assistance.
%
% Usage:
%   process_lws_folder('path/to/your/lws_folder'); % Regular processing
%   process_lws_folder('path/to/your/lws_folder', true); % With FFT smoothing
%   process_lws_folder('path/to/your/lws_folder', true, 8); % With custom smoothing factor

function process_lws_folder(inputDir, interpolate, interpolation_factor)

    if nargin < 1
        error('Usage: process_lws_folder(inputDir, [interpolate], [interpolation_factor])');
    end
    
    % Set defaults (always interpolate going forward)
    if nargin < 2 || isempty(interpolate)
        interpolate = true;
    end
    if nargin < 3 || isempty(interpolation_factor)
        interpolation_factor = 4;
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
    elseif strcmp(newLwsDir, inputDir)
        fprintf('Directory already named LWS, processing...\n');
    else
        fprintf('Directory %s already exists, moving files over\n', newLwsDir);
        movefile(fullfile(inputDir, '*'), newLwsDir);
        % Remove empty original directory if different from target
        if exist(inputDir, 'dir') && ~strcmp(inputDir, newLwsDir)
            rmdir(inputDir);
        end
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

    % Step 1: Convert .lws files to .csv (always interpolated)
    fprintf('\nStep 1: Converting .lws files to .csv format with FFT interpolation (factor: %d)...\n', interpolation_factor);
    fprintf('Input LWS directory: %s\n', inputDir);
    fprintf('Output CSV directory: %s\n', outputCsvDir);
    try
        lws_to_csv(inputDir, outputCsvDir, interpolation_factor); % unified
        fprintf('Successfully converted .lws files to interpolated .csv files.\n');
    catch ME
        fprintf('ERROR during LWS to CSV conversion: %s\n', ME.message);
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
        fprintf('Plot generation failed, but CSV conversion was successful.\n');
    end

    % Summary
    fprintf('\nProcessing complete for folder: %s (FFT interpolation factor: %d)\n', inputDir, interpolation_factor);

end 