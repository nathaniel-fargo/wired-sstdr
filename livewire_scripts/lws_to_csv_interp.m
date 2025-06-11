% lws_to_csv_interp.m
%
% Converts LWS files to CSV format with FFT-based interpolation applied to waveform data.
% The smoothing is performed before normalization to preserve signal characteristics.
%
% Author: Nathaniel Fargo
% Date: 2025-01-27
% Org: U of U WIRED 
%
% Developed partially with AI assistance.
%
% Usage:
%   lws_to_csv_interp('my_data.lws', 'my_data.csv');
%   lws_to_csv_interp('my_data.lws', 'my_data.csv', 4); % Custom interpolation factor
%   lws_to_csv_interp('my_data.lws', 'output_folder/');
%   lws_to_csv_interp('input_lws_directory/', 'output_csv_directory/');
%   lws_to_csv_interp('input_lws_directory/', 'output_csv_directory/', 8); % Higher smoothing

function lws_to_csv_interp(inputPath, outputPath, interpolation_factor)

    % Set default interpolation factor if not provided
    if nargin < 3
        interpolation_factor = 4; % Default interpolation factor
    end

    if isfolder(inputPath)
        % Input is a directory
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        files = dir(fullfile(inputPath, '*.lws'));
        if isempty(files)
            disp(['No .lws files found in directory: ' inputPath]);
            return;
        end
        disp(['Starting batch processing of ' num2str(length(files)) ' .lws files with smoothing (factor: ' num2str(interpolation_factor) ')...']);
        for i = 1:length(files)
            inputFile = fullfile(inputPath, files(i).name);
            [~, name, ~] = fileparts(files(i).name);
            outputFile = fullfile(outputPath, [name '_smooth.csv']);
            disp(['Processing ' inputFile ' -> ' outputFile]);
            try
                process_single_lws_smooth_local(inputFile, outputFile, interpolation_factor);
            catch ME
                disp(['Error processing file ' inputFile ': ' ME.message]);
                fprintf('Stack trace:\n%s\n', ME.getReport('extended', 'hyperlinks','off'));
            end
        end
        disp('Batch processing complete.');
    elseif isfile(inputPath)
        % Input is a single file
        [~, ~, ext] = fileparts(inputPath);
        if ~strcmpi(ext, '.lws')
            error('Input file must be a .lws file.');
        end
        
        [~, inName, ~] = fileparts(inputPath);

        % Determine output file path
        if endsWith(outputPath, '/') || endsWith(outputPath, '\\') || exist(outputPath, 'dir') % Check for trailing slash or if it's an existing directory
            % outputPath is a directory
            if ~exist(outputPath, 'dir')
                mkdir(outputPath);
            end
            outputFile = fullfile(outputPath, [inName '_smooth.csv']);
        else
            % outputPath is a filename
            [parentDir, outNameOnly, outExtOnly] = fileparts(outputPath);
            if ~isempty(parentDir) && ~exist(parentDir, 'dir')
                mkdir(parentDir);
            end
            if ~strcmpi(outExtOnly, '.csv')
                if isempty(outNameOnly) 
                     outputFile = fullfile(parentDir, [inName '_smooth.csv']);
                else
                     outputFile = fullfile(parentDir, [outNameOnly '_smooth.csv']);
                end
            else
                outputFile = outputPath; % Already path/name.csv
            end
        end
        
        disp(['Processing ' inputPath ' -> ' outputFile ' with smoothing (factor: ' num2str(interpolation_factor) ')']);
        try
            process_single_lws_smooth_local(inputPath, outputFile, interpolation_factor);
            disp('Processing complete.');
        catch ME
            disp(['Error processing file ' inputPath ': ' ME.message]);
            fprintf('Stack trace:\n%s\n', ME.getReport('extended', 'hyperlinks','off'));
            rethrow(ME);
        end
    else
        error(['Input path "' inputPath '" is not a valid file or directory.']);
    end
end

function process_single_lws_smooth_local(inputFile, outputFile, interpolation_factor)
    fid = fopen(inputFile, 'r');
    if fid == -1
        error(['Cannot open file: ' inputFile]);
    end
    
    % Read header lines
    line1_str = fgetl(fid);
    line2_str = fgetl(fid);
    line3_str = fgetl(fid);
    line4_str = fgetl(fid);
    
    if ~ischar(line1_str) || ~ischar(line2_str) || ~ischar(line3_str) || ~ischar(line4_str)
        fclose(fid);
        error(['File ' inputFile ' is too short or not in expected LWS format (missing header lines).']);
    end

    % Parse header lines
    header1_data = regexp(line1_str, 'SerialNumber=([^,]+), SoftwareVersion=([^,]+), Date=(.+)', 'tokens');
    serialNumber = ''; softwareVersion = ''; fileDate = '';
    if ~isempty(header1_data)
        serialNumber = strtrim(header1_data{1}{1});
        softwareVersion = strtrim(header1_data{1}{2});
        fileDate = strtrim(header1_data{1}{3});
    end
    
    header2_data = regexp(line2_str, 'Modulation=([^,]+), EndOffset=([^,]+), Resolution=([^,]+)', 'tokens');
    modulation = ''; endOffset = ''; resolution = '';
    if ~isempty(header2_data)
        modulation = strtrim(header2_data{1}{1});
        endOffset = strtrim(header2_data{1}{2});
        resolution = strtrim(header2_data{1}{3});
    end

    cableName = local_extractValue(line3_str, 'CableName');
    vop = local_extractValue(line3_str, 'VOP');
    preferredFrequency = local_extractValue(line3_str, 'PreferredFrequency');
    units = local_extractValue(line3_str, 'Units');
    minFaultDistance = local_extractValue(line3_str, 'MinFaultDistance');
    maxFaultDistance = local_extractValue(line3_str, 'MaxFaultDistance');

    selectedFrequency = local_extractValue(line4_str, 'SelectedFrequency');
    distanceAtSelection = local_extractValue(line4_str, 'Distance');

    csv_data_cell_array = {}; 
    
    csv_header = { ...
        'SerialNumber', 'SoftwareVersion', 'Date', 'Modulation', 'EndOffset', 'Resolution', ...
        'CableName', 'VOP', 'PreferredFrequency', 'Units', 'MinFaultDistance', 'MaxFaultDistance', ...
        'SelectedFrequencyAtAcquisition', 'DistanceAtAcquisition', ...
        'MeasurementFrequency', 'ZeroIndex', 'MinimumFaultOffset', 'MinimumFaultAmplitude', ...
        'FaultIndex', 'FaultAmplitude', 'UnitsPerSample', 'DataType', 'DataIndex', 'Value' ...
    };
    
    common_row_prefix = { ...
        serialNumber, softwareVersion, fileDate, modulation, endOffset, resolution, ...
        cableName, vop, preferredFrequency, units, minFaultDistance, maxFaultDistance, ...
        selectedFrequency, distanceAtSelection ...
    };

    file_line_number = 4;
    
    while true
        freq_line_str = fgetl(fid);
        file_line_number = file_line_number + 1;

        if ~ischar(freq_line_str) 
            break;
        end
        freq_line_str = strtrim(freq_line_str);
        if isempty(freq_line_str) 
            continue;
        end
        
        measFreq = local_extractValue(freq_line_str, 'Frequency');
        if isempty(measFreq)
            disp(['Warning: Could not parse "Frequency=" from line ' num2str(file_line_number) ' in ' inputFile ': "' freq_line_str '". Assuming end of data blocks or malformed file.']);
            break; 
        end
        
        zeroIdx = local_extractValue(freq_line_str, 'ZeroIndex');
        minFaultOffset = local_extractValue(freq_line_str, 'MinimumFaultOffset');
        minFaultAmp = local_extractValue(freq_line_str, 'MinimumFaultAmplitude');
        faultIdx = local_extractValue(freq_line_str, 'FaultIndex');
        faultAmp = local_extractValue(freq_line_str, 'FaultAmplitude');
        unitsPerSample = local_extractValue(freq_line_str, 'UnitsPerSample');
        
        current_freq_data_prefix = { ...
            measFreq, zeroIdx, minFaultOffset, minFaultAmp, ...
            faultIdx, faultAmp, unitsPerSample ...
        };
        
        % Process calibration data (unchanged - no smoothing for calibration)
        calib_line_str = fgetl(fid);
        file_line_number = file_line_number + 1;
        if ~ischar(calib_line_str) || ~startsWith(lower(strtrim(calib_line_str)), 'calibration:')
            disp(['Warning: Expected Calibration line after frequency data (line ' num2str(file_line_number-1) '), but got: "' calib_line_str ...
                  '" (or EOF) at input line ' num2str(file_line_number) ' in ' inputFile '. Stopping parse for this file.']);
            break; 
        end
        calib_data_str_parts = strsplit(strtrim(regexprep(calib_line_str, 'Calibration:\s*', '', 'ignorecase')), '\t');
        calib_data_num = str2double(calib_data_str_parts);
        for k = 1:length(calib_data_num)
            if ~isempty(calib_data_str_parts{k}) && ~isnan(calib_data_num(k)) 
                csv_data_cell_array(end+1, :) = [common_row_prefix, current_freq_data_prefix, {'Calibration'}, k, calib_data_num(k)];
            end
        end
        
        % Process waveform data with smoothing
        waveform_line_str = fgetl(fid);
        file_line_number = file_line_number + 1;
        if ~ischar(waveform_line_str) || ~startsWith(lower(strtrim(waveform_line_str)), 'waveform:')
            disp(['Warning: Expected Waveform line after calibration data (line ' num2str(file_line_number-1) '), but got: "' waveform_line_str ...
                  '" (or EOF) at input line ' num2str(file_line_number) ' in ' inputFile '. Stopping parse for this file.']);
            break;
        end
        waveform_data_str_parts = strsplit(strtrim(regexprep(waveform_line_str, 'Waveform:\s*', '', 'ignorecase')), '\t');
        waveform_data_num = str2double(waveform_data_str_parts);
        
        % Filter out NaN values and get valid data
        valid_indices = ~isnan(waveform_data_num) & ~cellfun(@isempty, waveform_data_str_parts);
        if sum(valid_indices) == 0
            disp(['Warning: No valid waveform data found for frequency ' measFreq ' in ' inputFile]);
            continue;
        end
        
        original_data = waveform_data_num(valid_indices);
        original_indices = find(valid_indices);
        
        % Apply FFT-based smoothing to waveform data
        try
            smoothed_data = apply_fft_smoothing(original_data, interpolation_factor);
            
            % Generate new indices for smoothed data
            original_spacing = 1; % Assume unit spacing in original indices
            new_spacing = original_spacing / interpolation_factor;
            smoothed_indices = (0:length(smoothed_data)-1) * new_spacing + original_indices(1);
            
            % Add smoothed waveform data to CSV
            for k = 1:length(smoothed_data)
                csv_data_cell_array(end+1, :) = [common_row_prefix, current_freq_data_prefix, {'Waveform'}, smoothed_indices(k), smoothed_data(k)];
            end
            
            % disp(['Applied smoothing to frequency ' measFreq ': ' num2str(length(original_data)) ' -> ' num2str(length(smoothed_data)) ' points']);
            
        catch ME_smooth
            disp(['Warning: Smoothing failed for frequency ' measFreq ' in ' inputFile ': ' ME_smooth.message]);
            disp('Using original data without smoothing.');
            
            % Fall back to original data
            for k = 1:length(original_data)
                if ~isnan(original_data(k))
                    csv_data_cell_array(end+1, :) = [common_row_prefix, current_freq_data_prefix, {'Waveform'}, original_indices(k), original_data(k)];
                end
            end
        end
    end
    
    fclose(fid);

    
    
    if isempty(csv_data_cell_array)
        disp(['Warning: No data parsed from ' inputFile '. CSV file ' outputFile ' will be created with headers only.']);
        T = cell2table(cell(0,length(csv_header)), 'VariableNames', csv_header);
    else
        if size(csv_data_cell_array,2) == length(csv_header)
            % Ensure numeric columns are numeric before creating table for robustness
            % DataIndex is col 23, Value is col 24
            for r = 1:size(csv_data_cell_array,1)
                if ischar(csv_data_cell_array{r, end-1}) % DataIndex
                    csv_data_cell_array{r, end-1} = str2double(csv_data_cell_array{r, end-1});
                end
                if ischar(csv_data_cell_array{r, end})   % Value
                    csv_data_cell_array{r, end} = str2double(csv_data_cell_array{r, end});
                end
            end
            T = cell2table(csv_data_cell_array, 'VariableNames', csv_header);
        else
            error(['Mismatch between number of columns in parsed data (' num2str(size(csv_data_cell_array,2)) ...
                   ') and expected CSV header columns (' num2str(length(csv_header)) ') for file ' inputFile '.']);
        end
    end
    
    try
        writetable(T, outputFile);
        disp(['Smoothed CSV file written: ' outputFile]);
    catch ME_write
        error(['Failed to write CSV file: ' outputFile '. Error: ' ME_write.message]);
    end
end

function smoothed_data = apply_fft_smoothing(data, interpolation_factor)
    % Apply FFT-based interpolation with zero padding (from plot_smooth_livewire.m)
    
    % Convert to row vector if needed
    if size(data, 1) > size(data, 2)
        data = data';
    end
    
    N_original = length(data);
    
    % Need at least 4 points for meaningful FFT interpolation
    if N_original < 4
        warning('Not enough data points for FFT interpolation. Using original data.');
        smoothed_data = data;
        return;
    end
    
    % Apply FFT-based interpolation with zero padding
    Y_fft = fft(data);
    N_new = N_original * interpolation_factor;
    
    % Zero pad in frequency domain
    Y_padded = zeros(1, N_new);
    half_N = floor(N_original/2);
    
    % Copy positive frequencies
    Y_padded(1:half_N+1) = Y_fft(1:half_N+1);
    
    % Copy negative frequencies
    if mod(N_original, 2) == 0
        % Even length: handle Nyquist frequency
        Y_padded(N_new-half_N+2:N_new) = Y_fft(half_N+2:N_original);
    else
        % Odd length
        Y_padded(N_new-half_N+1:N_new) = Y_fft(half_N+2:N_original);
    end
    
    % Scale to maintain energy
    Y_padded = Y_padded * interpolation_factor;
    
    % Inverse FFT to get interpolated signal
    smoothed_data = real(ifft(Y_padded));
end

function value = local_extractValue(line_str, key_name)
    value = ''; 
    
    if strcmpi(key_name, 'CableName')
        expr_start = [key_name '='];
        idx_start = regexpi(line_str, expr_start); % case insensitive search for key
        if ~isempty(idx_start)
            val_start_idx = idx_start(1) + length(expr_start);
            potential_next_delimiters_keys = {'VOP', 'PreferredFrequency', 'Units', 'MinFaultDistance', 'MaxFaultDistance'};
            end_val_idx = length(line_str);
            
            for i=1:length(potential_next_delimiters_keys)
                delim_pattern = [',\s*' potential_next_delimiters_keys{i} '=']; % comma, optional space, key=
                delim_match_info = regexpi(line_str, delim_pattern); % case insensitive
                
                if ~isempty(delim_match_info)
                    if delim_match_info(1) > val_start_idx
                        end_val_idx = min(end_val_idx, delim_match_info(1) -1); 
                    end
                end
            end
            value = strtrim(line_str(val_start_idx:end_val_idx));
        end
    else
        % Pattern: key_name=value_content (captures value_content until next comma or EOL)
        pattern = [key_name '\s*=\s*([^,]+?)(?:\s*,|$)']; % Non-greedy match, ends at comma or EOL
        tokens = regexp(line_str, pattern, 'tokens', 'once', 'ignorecase');
        if ~isempty(tokens)
            value = strtrim(tokens{1});
        else
            % Fallback for last parameter on the line (no trailing comma)
            pattern_end = [key_name '\s*=\s*(.*)'];
            tokens_end = regexp(line_str, pattern_end, 'tokens', 'once', 'ignorecase');
            if ~isempty(tokens_end)
                value = strtrim(tokens_end{1});
            end
        end
    end
end 