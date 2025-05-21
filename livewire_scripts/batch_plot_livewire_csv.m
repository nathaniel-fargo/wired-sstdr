% batch_plot_livewire_csv.m
%
% Processes all CSV files in a specified input directory, generates plots
% using plot_livewire_csv, and saves them as PNG files in an output directory.
%
% Author: Nathaniel Fargo
% Date: 2025-05-21
% Org: U of U WIRED
%
% Usage:
%   batch_plot_livewire_csvs('path/to/csv_folder', 'path/to/png_output_folder');


function batch_plot_livewire_csv(inputCsvDir, outputPngDir)

    if nargin < 2
        error('Usage: batch_plot_livewire_csv(inputCsvDir, outputPngDir)');
    end

    if ~exist(inputCsvDir, 'dir')
        error('Input CSV directory not found: %s', inputCsvDir);
    end

    if ~exist(outputPngDir, 'dir')
        fprintf('Output directory %s does not exist. Creating it...\\n', outputPngDir);
        mkdir(outputPngDir);
    end

    csvFiles = dir(fullfile(inputCsvDir, '*.csv'));
    if isempty(csvFiles)
        disp(['No CSV files found in ', inputCsvDir]);
        return;
    end

    fprintf('Found %d CSV files to process.\\n', length(csvFiles));

    for i = 1:length(csvFiles)
        currentCsvFile = csvFiles(i).name;
        fullCsvPath = fullfile(inputCsvDir, currentCsvFile);
        
        [~, baseName, ~] = fileparts(currentCsvFile);
        outputPngPath = fullfile(outputPngDir, [baseName, '.png']);
        
        fprintf('Processing %s -> %s\\n', fullCsvPath, outputPngPath);
        
        try
            % Call plot_livewire_csv with the output path.
            % We pass an empty cell {} for specificFrequencies to use default plotting logic.
            plot_livewire_csv(fullCsvPath, {}, outputPngPath);
            fprintf('Successfully generated plot: %s\\n', outputPngPath);
        catch ME
            fprintf('ERROR processing %s: %s\\n', fullCsvPath, ME.message);
            % Optionally, display more error details
            % disp(ME.getReport());
        end
    end

    disp('Batch processing complete.');

end 