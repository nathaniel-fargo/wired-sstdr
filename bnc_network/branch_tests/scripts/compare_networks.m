% compare_networks.m
%
% Compares waveform measurements for pairs of power-grid networks that differ by at most a user-specified number of wires (diffLevel).
% For each qualifying pair, selects the lowest mutual preferred measurement frequency, loads and normalizes the waveform, interpolates onto a common grid, and computes difference metrics.
% Optionally saves overlay plots and writes a summary CSV.
%
% Author: Nathaniel Fargo
% Date: 2025-05-25
% Org: U of U WIRED
%
% Usage:
%   compare_networks();
%   compare_networks(networkCsvPath, dataDirPath, diffLevel, generatePlots, fixedFrequency);

function compare_networks(networkCsvPath, dataDirPath, diffLevel, generatePlots, fixedFrequency)

%% ---------------------------- input parsing ---------------------------- %%
if nargin < 1 || isempty(networkCsvPath)
    networkCsvPath = "bnc_network/branch_tests/2025-05-22/networks_n.csv";
end
if nargin < 2 || isempty(dataDirPath)
    % Default to a LiveWire CSV directory that is *relative to the network CSV*
    dataDirPath = "LiveWire/Static/CSV";
end
if nargin < 3 || isempty(diffLevel)
    diffLevel = 1;
end
if nargin < 4 || isempty(generatePlots)
    generatePlots = false;
end
if nargin < 5
    fixedFrequency = "";
end

% If the path given (or defaulted) is *relative*, interpret it relative to
% the directory that contains the network list CSV, NOT the current MATLAB
% working directory.
if ~isAbsolutePath(dataDirPath)
    dataDirPath = fullfile(fileparts(networkCsvPath), dataDirPath);
end

assert(isfile(networkCsvPath),   "Network CSV file not found: %s", networkCsvPath);
assert(isfolder(dataDirPath),    "Data directory not found: %s",   dataDirPath);

%% ------------------------- read network catalogue ---------------------- %%
netTbl = readtable(networkCsvPath, "TextType", "string");
requiredCols = ["ID","Network"];
if ~all(ismember(requiredCols, netTbl.Properties.VariableNames))
    error("Network CSV must contain columns: %s", strjoin(requiredCols, ", "));
end

numNets = height(netTbl);

% Parse every network once and cache tokens / depth per network.
netTokens  = cell(numNets,1);   % cell{idx} -> string array of tokens
netDepths  = cell(numNets,1);   % matching depth values
for r = 1:numNets
    [toks, deps]    = parseNetworkString(netTbl.Network(r));
    netTokens{r}    = toks;
    netDepths{r}    = deps; %#ok<NASGU>  % kept in case depth reporting needed
end


%% ----------------------- prepare output containers --------------------- %%
resultsCell = cell(0, 9);  % Initialize with 0 rows, 9 columns

%% ----------------------- waveform cache (avoid IO) --------------------- %%
waveCache   = containers.Map();   % key = network ID -> struct with fields:
                                  %   .freqs         (string array)
                                  %   .preferredFreq (string)
                                  %   .data          (struct of frequency -> struct(x,y))

%% ---------------------- pair-wise comparisons -------------------------- %%
for i = 1:numNets-1
    idA   = netTbl.ID(i);
    toksA = netTokens{i};
    for j = i+1:numNets
        idB   = netTbl.ID(j);
        toksB = netTokens{j};

        % --- difference metric (# wires not shared) -------------------- %
        diffWires = setdiff(union(toksA, toksB), intersect(toksA, toksB));
        nDiff     = numel(diffWires);
        
        % --- compute minimum depth of differing wires --- %
        depthsA = netDepths{i};
        depthsB = netDepths{j};
        minDepths = nan(1, nDiff);
        for d = 1:nDiff
            w = diffWires(d);
            idxA = find(toksA == w, 1);
            idxB = find(toksB == w, 1);
            depthVals = [];
            if ~isempty(idxA), depthVals(end+1) = depthsA(idxA); end
            if ~isempty(idxB), depthVals(end+1) = depthsB(idxB); end
            if ~isempty(depthVals)
                minDepths(d) = min(depthVals);
            end
        end
        if isempty(minDepths(~isnan(minDepths)))
            minDiffDepth = NaN;
        else
            minDiffDepth = min(minDepths(~isnan(minDepths)));
        end

        if nDiff > diffLevel
            continue;  % skip – too different
        end

        try
            % ---------------- load waveform info ---------------------- %
            dataA = loadNetwork(idA, netTbl.Network(i), dataDirPath, waveCache);
            dataB = loadNetwork(idB, netTbl.Network(j), dataDirPath, waveCache);
        catch ME
            warning("Failed loading data for pair %s/%s: %s", idA, idB, ME.message);
            continue;
        end

        % ------------- choose frequency (fixed or preferred) ------------- %
        if fixedFrequency ~= ""
            if any(dataA.freqs == fixedFrequency) && any(dataB.freqs == fixedFrequency)
                freqChosen = fixedFrequency;
            else
                warning("Fixed frequency '%s' not available for pair %s/%s - skipping", fixedFrequency, idA, idB);
                continue;
            end
        else
            % ------------- choose common frequency (lowest preferred) ------------- %
            freqsCommon = intersect(dataA.freqs, dataB.freqs);
            if isempty(freqsCommon)
                warning("No common measurement frequency for %s & %s - skipping", idA, idB);
                continue;
            end

            prefA = dataA.preferredFreq;
            prefB = dataB.preferredFreq;
            inA = any(freqsCommon == prefA);
            inB = any(freqsCommon == prefB);
            chosen = "";
            if inA && inB
                % Both present: choose the lower (in Hz)
                hzA = parseFreqStr(prefA);
                hzB = parseFreqStr(prefB);
                if hzA <= hzB
                    chosen = prefA;
                else
                    chosen = prefB;
                end
            elseif inA
                chosen = prefA;
            elseif inB
                chosen = prefB;
            end
            if chosen ~= ""
                freqChosen = chosen;
            else
                % Fallback: choose the lowest mutual frequency
                freqHzVals = arrayfun(@parseFreqStr, cellstr(freqsCommon));
                [~, idxMin] = min(freqHzVals);
                freqChosen  = freqsCommon(idxMin);
            end
        end

        wfA = dataA.data(freqChosen);
        wfB = dataB.data(freqChosen);

        % Ensure we have non-empty waveforms
        if isempty(wfA.x) || isempty(wfB.x)
            warning("Empty waveform at %s for %s or %s", freqChosen, idA, idB);
            continue;
        end

        % ------------- interpolate onto common grid -------------------- %
        xCommon = unique([wfA.x; wfB.x]);
        yA      = interp1(wfA.x, wfA.y, xCommon, "linear", "extrap");
        yB      = interp1(wfB.x, wfB.y, xCommon, "linear", "extrap");
        yDiff   = yA - yB;

        % -------------- metrics --------------------------------------- %
        peakMag       = max(abs(yDiff));
        areaUnsquared = trapz(xCommon, yDiff);
        areaSquared   = trapz(xCommon, yDiff.^2);

        % -------------- store result ---------------------------------- %
        resultsCell(end+1, :) = {idA, idB, strjoin(diffWires, ";"), nDiff, ...
                                 freqChosen, peakMag, areaUnsquared, areaSquared, minDiffDepth};

        % -------------- optional plot --------------------------------- %
        if generatePlots
            saveDifferencePlot(idA, idB, xCommon, yA, yB, yDiff, freqChosen, dataDirPath);
        end
    end
end

%% --------------------------- write CSV -------------------------------- %%
if ~isempty(resultsCell)
    outTbl = cell2table(resultsCell, 'VariableNames', ...
        {'NetworkA','NetworkB','DifferingWires','NumWireDiff', ...
         'Frequency','PeakMag','AreaUnsquared','AreaSquared','MinDiffDepth'});

    [outDir, baseName, ~] = fileparts(networkCsvPath);
    outFile = fullfile(outDir, baseName + "_analysis.csv");
    writetable(outTbl, outFile);
    fprintf("Saved pair-analysis results to %s\n", outFile);
else
    fprintf("No network pairs met the diff-level criterion (%d).\n", diffLevel);
end

end  % main function

%% ====================================================================== %%
function [tokens, depths] = parseNetworkString(netStr)
%PARSE_NETWORK_STRING  Extract wire identifiers and their depths from a
% network definition string such as "{D03{E00[O]}...}".
% The grammar (see README) is a series of nested braces "{<wire> ...}".
% We simply scan the string character-by-character, maintaining a depth
% counter and recording any occurrence of an identifier matching the regex
% /[A-Z][0-9]{2}/ that is *outside* square-bracket termination annotations.

    if ismissing(netStr) || strlength(netStr)==0
        tokens = string.empty;
        depths = double.empty;
        return;
    end

    s      = char(netStr);      % work with char vector for speed
    n      = length(s);
    d      = -1;                % depth before seeing first '{' will become 0
    tokens = string.empty;
    depths = double.empty;
    inTerm = false;             % inside [...] termination annotation

    i = 1;
    while i <= n
        ch = s(i);
        switch ch
            case '{'
                d = d + 1;
                i = i + 1;
            case '}'
                d = max(d-1,0);
                i = i + 1;
            case '['
                inTerm = true;
                i = i + 1;
            case ']'
                inTerm = false;
                i = i + 1;
            otherwise
                if ~inTerm && ch>='A' && ch<='Z' && i+2 <= n && isstrprop(s(i+1),'digit') && isstrprop(s(i+2),'digit')
                    tok = s(i:i+2);
                    tokens(end+1) = string(tok); %#ok<AGROW>
                    depths(end+1) = d;          %#ok<AGROW>
                    i = i + 3;  % skip past token
                else
                    i = i + 1;  % move on
                end
        end
    end
end

%% ====================================================================== %%
function dataStruct = loadNetwork(netID, netString, dataDir, cache)
% Retrieve (and cache) waveform data for a network.  The cache avoids
% re-reading the same CSV multiple times.

    if isKey(cache, netID)
        dataStruct = cache(netID);
        return;
    end

    % Guess measurement CSV filename – assume first file that starts with ID
    fPattern = fullfile(dataDir, netID + "*.csv");
    fList    = dir(fPattern);
    if isempty(fList)
        error("No measurement file matching %s", fPattern);
    end
    csvFile  = fullfile(fList(1).folder, fList(1).name);

    % Read once and preprocess
    opts = detectImportOptions(csvFile);
    % Ensure critical columns are read as strings, similar to plot_livewire_csv.m
    stringColumns = {"MeasurementFrequency", "PreferredFrequency", "DataType", "Units", ...
                     "UnitsPerSample", "ZeroIndex", "DistanceAtAcquisition"};
    for k = 1:length(stringColumns)
        colName = stringColumns{k};
        if any(strcmp(opts.VariableNames, colName))
            opts = setvartype(opts, colName, 'string');
        end
    end
    dataTbl = readtable(csvFile, opts);

    % Waveform rows only
    wfRows  = strcmpi(dataTbl.DataType, "Waveform");
    wfTbl   = dataTbl(wfRows, :);

    if isempty(wfTbl)
        error("No waveform rows in %s", csvFile);
    end

    uniqueFreqs = unique(wfTbl.MeasurementFrequency, "stable");

    % Determine preferred frequency (use PreferredFrequency, else SelectedFrequencyAtAcquisition, else first available)
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
        f            = uniqueFreqs(k);
        freqRows     = wfTbl.MeasurementFrequency == f;
        [x,y]        = waveformToXY(wfTbl(freqRows,:));
        freqMap(char(f)) = struct('x',x,'y',y);
    end

    dataStruct = struct('freqs', uniqueFreqs, ...
                        'preferredFreq', prefFreq, ...
                        'data', freqMap);
    cache(netID) = dataStruct;  % store


end

%% ====================================================================== %%
function [x_ft, y_norm] = waveformToXY(tbl)
% Convert a table slice for one frequency into X (ft) and normalised Y.

    METERS_TO_FEET = 3.28084;

    % Y normalisation
    y_raw   = tbl.Value;
    maxAbs  = max(abs(y_raw));
    if maxAbs == 0
        y_norm = zeros(size(y_raw));
    else
        y_norm = y_raw / maxAbs;
    end

    % Distance in original units
    dataIdx      = tbl.DataIndex;
    upsRaw       = tbl.UnitsPerSample(1);
    if isstring(upsRaw) || ischar(upsRaw)
        ups = str2double(upsRaw);
    else
        ups = upsRaw;
    end
    if isnan(ups) || ups <= 0, ups = 1.0; end
    distOrig     = dataIdx .* ups; % same units as UPS

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

    zeroOffset   = zIdx * ups;  % same units as distOrig

    if isMetric
        distFt = (distOrig - zeroOffset) * METERS_TO_FEET;
    else
        distFt = (distOrig - zeroOffset);
    end

    x_ft = distFt;
end

%% ====================================================================== %%
function saveDifferencePlot(idA, idB, x, yA, yB, yDiff, freqStr, dataDir)
% Create and save a PNG that overlays yA, yB, and their difference.

    parentDir   = fileparts(dataDir); % one level above CSV
    diffDir     = fullfile(parentDir, "Plots", "Differences");
    if ~isfolder(diffDir), mkdir(diffDir); end

    fig = figure('Visible','off');
    plot(x, yA,  '-', 'DisplayName', idA);
    hold on;
    plot(x, yB,  '-', 'DisplayName', idB);
    plot(x, yDiff, '--', 'DisplayName', 'Difference');
    hold off;
    legend('Location', 'best');
    xlabel('Distance (ft)'); ylabel('Normalised Value');
    title(sprintf('%s vs %s – %s', idA, idB, freqStr));
    grid on;

    fname = sprintf('%s_%s_diff.png', idA, idB);
    saveas(fig, fullfile(diffDir, fname));
    close(fig);
end

%% ---------------------------------------------------------------------- %%
function tf = isAbsolutePath(p)
% Return true if p is an absolute path on this platform.
    if isempty(p); tf = false; return; end
    p = char(p);
    if ispc
        tf = length(p) >= 2 && p(2) == ':';           % e.g. C:\...
    else
        tf = p(1) == '/';                             % POSIX absolute path
    end
end

% Place this helper at the end of the file as a subfunction
function hz = parseFreqStr(freqStr)
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
