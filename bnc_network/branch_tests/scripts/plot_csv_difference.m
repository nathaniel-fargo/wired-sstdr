% plot_csv_difference.m
%
% Standalone function to plot the difference between two LiveWire CSV measurements.
% Takes two CSV files, loads their waveform data, normalizes and interpolates them onto 
% a common grid, then creates an overlay plot showing both waveforms and their difference.
%
% Author: Nathaniel Fargo
% Date: 2025-01-07
% Org: U of U WIRED
%
% Usage:
%   plot_csv_difference(csvFileA, csvFileB);
%   plot_csv_difference(csvFileA, csvFileB, frequency);
%   plot_csv_difference(csvFileA, csvFileB, frequency, outputDir);
%   plot_csv_difference(csvFileA, csvFileB, frequency, outputDir, nameA, nameB);

function plot_csv_difference(csvFileA, csvFileB, frequency, outputDir, nameA, nameB)

%% ---------------------------- input parsing ---------------------------- %%
if nargin < 2
    error("Two CSV file paths are required");
end
if nargin < 3 || isempty(frequency)
    frequency = ""; % Auto-select preferred frequency
end
if nargin < 4 || isempty(outputDir)
    outputDir = pwd; % Current directory
end
if nargin < 5 || isempty(nameA)
    [~, nameA, ~] = fileparts(csvFileA);
end
if nargin < 6 || isempty(nameB)
    [~, nameB, ~] = fileparts(csvFileB);
end

% Validate input files
assert(isfile(csvFileA), "CSV file A not found: %s", csvFileA);
assert(isfile(csvFileB), "CSV file B not found: %s", csvFileB);

% Ensure output directory exists
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% --------------------------- load CSV data ----------------------------- %%
fprintf("Loading %s...\n", csvFileA);
dataA = loadCsvData(csvFileA);
fprintf("Loading %s...\n", csvFileB);
dataB = loadCsvData(csvFileB);

%% ------------------------- choose frequency ---------------------------- %%
if frequency ~= ""
    % Use specified frequency
    if any(dataA.freqs == frequency) && any(dataB.freqs == frequency)
        freqChosen = frequency;
    else
        error("Frequency '%s' not available in both files", frequency);
    end
else
    % Auto-select frequency
    freqsCommon = intersect(dataA.freqs, dataB.freqs);
    if isempty(freqsCommon)
        error("No common measurement frequency between the two files");
    end
    
    % Try to use preferred frequencies
    prefA = dataA.preferredFreq;
    prefB = dataB.preferredFreq;
    inA = any(freqsCommon == prefA);
    inB = any(freqsCommon == prefB);
    
    if inA && inB
        % Both present: choose the lower (in Hz)
        hzA = parseFreqStr(prefA);
        hzB = parseFreqStr(prefB);
        if hzA <= hzB
            freqChosen = prefA;
        else
            freqChosen = prefB;
        end
    elseif inA
        freqChosen = prefA;
    elseif inB
        freqChosen = prefB;
    else
        % Fallback: choose the lowest mutual frequency
        freqHzVals = arrayfun(@parseFreqStr, cellstr(freqsCommon));
        [~, idxMin] = min(freqHzVals);
        freqChosen = freqsCommon(idxMin);
    end
end

fprintf("Using frequency: %s\n", freqChosen);

%% ------------------------- extract waveforms --------------------------- %%
wfA = dataA.data(freqChosen);
wfB = dataB.data(freqChosen);

% Ensure we have non-empty waveforms
if isempty(wfA.x) || isempty(wfB.x)
    error("Empty waveform at %s for one or both files", freqChosen);
end

%% ------------------- interpolate onto common grid ---------------------- %%
xCommon = unique([wfA.x; wfB.x]);
yA = interp1(wfA.x, wfA.y, xCommon, "linear", "extrap");
yB = interp1(wfB.x, wfB.y, xCommon, "linear", "extrap");
yDiff = yA - yB;

%% -------------------------- compute metrics ---------------------------- %%
peakMag = max(abs(yDiff));
areaNorm = trapz(xCommon, abs(yDiff));
areaSquared = trapz(xCommon, yDiff.^2);

fprintf("Difference Metrics:\n");
fprintf("  Peak Magnitude: %.4g\n", peakMag);
fprintf("  Area (L1 norm): %.4g\n", areaNorm);
fprintf("  Area Squared (L2^2): %.4g\n", areaSquared);

%% --------------------------- create plot ------------------------------- %%
fig = figure('Position', [100, 100, 800, 600]);

% Plot the waveforms and their difference
plot(xCommon, yA, 'Color', 'blue', 'LineWidth', 1, 'DisplayName', nameA); % Orange
hold on;
plot(xCommon, yB, 'Color', [0, 0.8, 0], 'LineWidth', 1, 'DisplayName', nameB); % Green
plot(xCommon, yDiff, '--', 'Color', [0.8, 0.557,  0], 'LineWidth', 1, 'DisplayName', 'Difference'); % Blue
hold off;

% Formatting
legend('Location', 'best');
xlabel('Distance (ft)');
ylabel('Normalized Value');
title(sprintf('%s vs %s at %s', nameA, nameB, freqChosen));
grid on;

% Add metrics as text annotation
metricsText = sprintf('Peak Mag: %.3g\nArea Norm: %.3g\nArea²: %.3g', ...
                     peakMag, areaNorm, areaSquared);
annotation('textbox', [0.02, 0.98, 0.25, 0.1], 'String', metricsText, ...
          'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'white', ...
          'EdgeColor', 'black');

%% --------------------------- display plot ------------------------------ %%
% Display the plot interactively for manual editing/cropping
figure(fig); % Bring to front and display
fprintf("Plot displayed. You can now crop, edit, or save manually.\n");
fprintf("Suggested filename: %s_vs_%s_diff_%s.png\n", nameA, nameB, strrep(freqChosen, ' ', '_'));

end

%% ====================================================================== %%
function dataStruct = loadCsvData(csvFile)
% Load and preprocess waveform data from a LiveWire CSV file

    % Read CSV with proper options
    opts = detectImportOptions(csvFile);
    
    % Ensure critical columns are read as strings
    stringColumns = {"MeasurementFrequency", "PreferredFrequency", "DataType", "Units", ...
                     "UnitsPerSample", "ZeroIndex", "DistanceAtAcquisition"};
    for k = 1:length(stringColumns)
        colName = stringColumns{k};
        if any(strcmp(opts.VariableNames, colName))
            opts = setvartype(opts, colName, 'string');
        end
    end
    
    dataTbl = readtable(csvFile, opts);
    
    % Filter to waveform rows only
    wfRows = strcmpi(dataTbl.DataType, "Waveform");
    wfTbl = dataTbl(wfRows, :);
    
    if isempty(wfTbl)
        error("No waveform rows found in %s", csvFile);
    end
    
    % Get unique frequencies
    uniqueFreqs = unique(wfTbl.MeasurementFrequency, "stable");
    
    % Determine preferred frequency
    prefFreq = "";
    if any(strcmpi(dataTbl.Properties.VariableNames, "PreferredFrequency"))
        prefCol = dataTbl.PreferredFrequency;
        idx = find(~ismissing(prefCol) & prefCol ~= "None", 1, 'first');
        if ~isempty(idx)
            prefFreq = prefCol(idx);
        end
    end
    if (prefFreq == "")
        if any(strcmpi(dataTbl.Properties.VariableNames, "SelectedFrequencyAtAcquisition"))
            selCol = dataTbl.SelectedFrequencyAtAcquisition;
            idx = find(~ismissing(selCol) & selCol ~= "None", 1, 'first');
            if ~isempty(idx)
                prefFreq = selCol(idx);
            end
        end
    end
    if (prefFreq == "")
        prefFreq = uniqueFreqs(1);
    end
    
    % Build per-frequency waveform structs
    freqMap = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(uniqueFreqs)
        f = uniqueFreqs(k);
        freqRows = wfTbl.MeasurementFrequency == f;
        [x, y] = waveformToXY(wfTbl(freqRows, :));
        freqMap(char(f)) = struct('x', x, 'y', y);
    end
    
    dataStruct = struct('freqs', uniqueFreqs, ...
                        'preferredFreq', prefFreq, ...
                        'data', freqMap);
end

%% ====================================================================== %%
function [x_ft, y_norm] = waveformToXY(tbl)
% Convert a table slice for one frequency into X (ft) and normalized Y.

    METERS_TO_FEET = 3.28084;
    
    % Y normalization
    y_raw = tbl.Value;
    maxAbs = max(abs(y_raw));
    if maxAbs == 0
        y_norm = zeros(size(y_raw));
    else
        y_norm = y_raw / maxAbs;
    end
    
    % Distance in original units
    dataIdx = tbl.DataIndex;
    upsRaw = tbl.UnitsPerSample(1);
    if isstring(upsRaw) || ischar(upsRaw)
        ups = str2double(upsRaw);
    else
        ups = upsRaw;
    end
    if isnan(ups) || ups <= 0, ups = 1.0; end
    distOrig = dataIdx .* ups; % same units as UPS
    
    % Units (Metric / Standard)
    unitsStr = tbl.Units(1);
    if ismissing(unitsStr) || unitsStr == "", unitsStr = "Standard"; end
    isMetric = strcmpi(unitsStr, "Metric");
    
    % Zero index offset
    zRaw = tbl.ZeroIndex(1);
    if isstring(zRaw) || ischar(zRaw)
        zIdx = str2double(zRaw);
    else
        zIdx = zRaw;
    end
    if isnan(zIdx), zIdx = 0; end
    
    zeroOffset = zIdx * ups; % same units as distOrig
    
    if isMetric
        distFt = (distOrig - zeroOffset) * METERS_TO_FEET;
    else
        distFt = (distOrig - zeroOffset);
    end
    
    x_ft = distFt;
end

%% ====================================================================== %%
function hz = parseFreqStr(freqStr)
% Parse frequency string to Hz value

    if iscell(freqStr)
        freqStr = freqStr{1};
    end
    freqStr = string(freqStr);
    freqStr = strtrim(freqStr);
    tokens = regexp(freqStr, '([\d.]+)\s*([kKmMgG]?[hH][zZ])', 'tokens', 'once');
    if isempty(tokens) || numel(tokens) < 2
        hz = NaN;
        warning('Could not parse frequency string: "%s"', freqStr);
        return;
    end
    val = str2double(tokens{1});
    unit = lower(tokens{2});
    switch unit
        case 'hz', hz = val;
        case 'khz', hz = val * 1e3;
        case 'mhz', hz = val * 1e6;
        case 'ghz', hz = val * 1e9;
        otherwise, hz = NaN;
    end
end 