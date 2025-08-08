% plot_data_with_baseline_and_lengths.m
% 
% Plots PeakMag vs. MinDiffDepth from a "trend" CSV, fits an exponential curve,
% plots an error level based on PeakMag from an "error" CSV,
% calculates and annotates their intersection,
% and shows a histogram of cable lengths from cable_logs.csv.
%
% Author: Nathaniel Fargo
% Date: 2025-06-02
% Org: U of U WIRED
%
% Usage:
%   plot_data_with_baseline_and_lengths('path/to/trend.csv', 'path/to/error.csv');
%   plot_data_with_baseline_and_lengths(); % prompts for files

function plot_data_with_baseline_and_lengths(trendCsvPath, errorCsvPath)
% DEPRECATED: use plot_data_with_baseline.m for the core analysis.

% --- Handle Input Arguments ---
if nargin < 1 || isempty(trendCsvPath)
    [f_trend,p_trend] = uigetfile('*.csv', 'Select TREND analysis CSV');
    if isequal(f_trend,0), disp('No trend file selected.'); return; end
    trendCsvPath = fullfile(p_trend,f_trend);
end

if nargin < 2 || isempty(errorCsvPath)
    [f_error,p_error] = uigetfile('*.csv', 'Select ERROR analysis CSV');
    if isequal(f_error,0), disp('No error file selected.'); return; end
    errorCsvPath = fullfile(p_error,f_error);
end

assert(isfile(trendCsvPath), 'Trend CSV file not found: %s', trendCsvPath);
assert(isfile(errorCsvPath), 'Error CSV file not found: %s', errorCsvPath);

% --- Load and Process Trend Data ---
trendTbl = readtable(trendCsvPath);
% Remove NaNs in MinDiffDepth and PeakMag
mask = ~isnan(trendTbl.MinDiffDepth) & ~isnan(trendTbl.PeakMag);
trendTbl = trendTbl(mask,:);
x_trend = trendTbl.MinDiffDepth;
y_trend_peakmag = trendTbl.PeakMag;

% --- Load and Process Error Data ---
errorTbl = readtable(errorCsvPath);
maskE = ~isnan(errorTbl.PeakMag);
errorTbl = errorTbl(maskE,:);
y_error_peakmag_values = errorTbl.PeakMag;
errorLevelPeakMag = mean(y_error_peakmag_values);

% --- Set up figure with two panels ---
fig = figure('Name','PeakMag Trend & Cable Lengths','Color','w','Position',[100 100 900 800]);
tiledlayout(2,1);

% --- Panel 1: Trend & Error Plot ---
nexttile;
hold on;
% Scatter trend data
scatter(x_trend, y_trend_peakmag, 36, 'filled', 'DisplayName', 'Trend Data');

% Exponential fit: y = a * exp(b*x)
fitMask = y_trend_peakmag > 0;
xFit = x_trend(fitMask);
yFit = y_trend_peakmag(fitMask);
if numel(xFit) >= 2
    p = polyfit(xFit, log(yFit), 1);
    a = exp(p(2)); b = p(1);
    xFitLine = linspace(min(x_trend), max(x_trend), 200);
    yExpFitLine = a * exp(b * xFitLine);
    plot(xFitLine, yExpFitLine, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Exp fit (y=%.2ge^{%.2gx})', a, b));
    % Intersection
    if a > 0 && errorLevelPeakMag > 0 && b ~= 0
        x_int = log(errorLevelPeakMag / a) / b;
        y_int = errorLevelPeakMag;
        plot(x_int, y_int, 'ko', 'MarkerSize',10,'MarkerFaceColor','k','DisplayName','Intersection');
    end
end
% Plot error level
x_limits = xlim;
if isempty(x_limits) || any(isnan(x_limits))
    x_limits = [min(x_trend)-0.5, max(x_trend)+0.5];
end
plot(x_limits, [errorLevelPeakMag, errorLevelPeakMag], 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Error Level = %.2g', errorLevelPeakMag));
xlim(x_limits);

xlabel('MinDiffDepth');
ylabel('PeakMag');
title('PeakMag vs. Branch Depth with Error Level');
grid on;
legend('show','Location','best');
hold off;

% --- Panel 2: Cable Length Histogram ---
nexttile;
% Locate cable_logs.csv relative to this script
thisFile = mfilename('fullpath');
scriptDir = fileparts(thisFile);
branchTestsDir = fileparts(scriptDir);
projectRoot = fileparts(branchTestsDir);
logsPath = fullfile(projectRoot, 'cable_measurements', 'cable_logs.csv');
assert(isfile(logsPath), 'Cable logs file not found: %s', logsPath);

logsTbl = readtable(logsPath);
% Parse lengths into numeric feet
lenFeet = nan(height(logsTbl),1);
for k = 1:height(logsTbl)
    lenFeet(k) = parseLengthFeet(logsTbl.Length{k});
end
histogram(lenFeet);
xlabel('Cable Length (ft)');
ylabel('Count');
title('Distribution of Cable Lengths');
grid on;

% --- Save Figure ---
[~,trendBase,~] = fileparts(trendCsvPath);
[~,errorBase,~] = fileparts(errorCsvPath);
outDir = fileparts(trendCsvPath);
if isempty(outDir), outDir = pwd; end
outFile = fullfile(outDir, sprintf('%s_vs_%s_PeakMag_with_lengths.png', trendBase, errorBase));
try
    saveas(fig, outFile);
    fprintf('Saved plot to %s\n', outFile);
catch ME
    fprintf('Could not save plot to %s: %s\n', outFile, ME.message);
    altFile = fullfile(pwd, sprintf('%s_vs_%s_PeakMag_with_lengths.png', trendBase, errorBase));
    try
        saveas(fig, altFile);
        fprintf('Saved plot to %s\n', altFile);
    catch ME2
        fprintf('Failed to save plot: %s\n', ME2.message);
    end
end

end  % main function

function ft = parseLengthFeet(lenStr)
% Converts a length string like '3ft 11in' to numeric feet
    ft = NaN;
    if ismissing(lenStr) || isempty(lenStr)
        return;
    end
    if iscell(lenStr), lenStr = lenStr{1}; end
    lenStr = char(lenStr);
    tokens = regexp(lenStr, '(?<ft>\d+)\s*ft(?:\s*(?<in>\d+)\s*in)?', 'names');
    if isempty(tokens)
        warning('Could not parse length: %s', lenStr);
        return;
    end
    ft = str2double(tokens.ft);
    if isfield(tokens,'in') && ~isempty(tokens.in)
        ft = ft + str2double(tokens.in)/12;
    end
end 