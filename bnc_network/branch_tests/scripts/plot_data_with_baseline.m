% plot_data_with_baseline.m
%
% Plots specified metric vs. MinDiffDepth from a "trend" CSV, fits an exponential curve,
% plots an error level based on the same metric from an "error" CSV,
% and calculates and annotates their intersection.
%
% Author: Nathaniel Fargo
% Date: 2025-05-28
% Org: U of U WIRED
%
% Usage:
%   plot_data_with_baseline('path/to/trend_analysis.csv', 'path/to/error_analysis.csv', 'PeakMag');
%   plot_data_with_baseline('path/to/trend_analysis.csv', 'path/to/error_analysis.csv', 'AreaNorm');
%   plot_data_with_baseline('path/to/trend_analysis.csv', 'path/to/error_analysis.csv', 'AreaSquared');
%   plot_data_with_baseline(); % prompts for files and metric

function plot_data_with_baseline(trendCsvPath, errorCsvPath, metric)

% --- Handle Input Arguments ---
if nargin < 1 || isempty(trendCsvPath)
    [f_trend,p_trend] = uigetfile('*.csv', 'Select TREND analysis CSV (e.g., networks_n_m_analysis.csv)');
    if isequal(f_trend,0), disp('No trend file selected.'); return; end
    trendCsvPath = fullfile(p_trend,f_trend);
end

if nargin < 2 || isempty(errorCsvPath)
    [f_error,p_error] = uigetfile('*.csv', 'Select ERROR analysis CSV (e.g., networks_v_analysis.csv)');
    if isequal(f_error,0), disp('No error file selected.'); return; end
    errorCsvPath = fullfile(p_error,f_error);
end

if nargin < 3 || isempty(metric)
    % Prompt user to select metric
    validMetrics = {'PeakMag', 'AreaNorm', 'AreaSquared'};
    [selection, ok] = listdlg('PromptString', 'Select metric to analyze:', ...
                              'SelectionMode', 'single', ...
                              'ListString', validMetrics);
    if ~ok, disp('No metric selected.'); return; end
    metric = validMetrics{selection};
end

% Validate metric
validMetrics = {'PeakMag', 'AreaNorm', 'AreaSquared'};
if ~ismember(metric, validMetrics)
    error('Invalid metric "%s". Must be one of: %s', metric, strjoin(validMetrics, ', '));
end

assert(isfile(trendCsvPath), 'Trend CSV file not found: %s', trendCsvPath);
assert(isfile(errorCsvPath), 'Error CSV file not found: %s', errorCsvPath);

% --- Load and Process Trend Data ---
trendTbl = readtable(trendCsvPath);
fprintf('Loaded trend data from: %s\n', trendCsvPath);
fprintf('Using metric: %s\n', metric);

% Remove rows with NaN in relevant columns for trend data
trendRelevantCols = {'MinDiffDepth', metric};
trendRowMask = true(height(trendTbl),1);
for k = 1:numel(trendRelevantCols)
    if ismember(trendRelevantCols{k}, trendTbl.Properties.VariableNames)
        trendRowMask = trendRowMask & ~isnan(trendTbl.(trendRelevantCols{k}));
    else
        error('Column "%s" not found in trend CSV: %s', trendRelevantCols{k}, trendCsvPath);
    end
end
trendTbl = trendTbl(trendRowMask, :);

if isempty(trendTbl)
    error('No valid data rows found in trend CSV after NaN removal for MinDiffDepth and %s.', metric);
end

x_trend = trendTbl.MinDiffDepth;
y_trend_metric = trendTbl.(metric);

fprintf('\nSummary statistics for Trend %s:\n', metric);
fprintf('  Mean:   %.4g\n', mean(y_trend_metric, 'omitnan'));
fprintf('  Median: %.4g\n', median(y_trend_metric, 'omitnan'));
fprintf('  Std:    %.4g\n', std(y_trend_metric, 'omitnan'));

% --- Load and Process Error Data ---
errorTbl = readtable(errorCsvPath);
fprintf('Loaded error data from: %s\n', errorCsvPath);

% Remove rows with NaN in the specified metric for error data
if ~ismember(metric, errorTbl.Properties.VariableNames)
    error('Column "%s" not found in error CSV: %s', metric, errorCsvPath);
end
errorRowMask = ~isnan(errorTbl.(metric));
errorTbl = errorTbl(errorRowMask, :);

if isempty(errorTbl)
    error('No valid data rows found in error CSV after NaN removal for %s.', metric);
end

y_error_metric_values = errorTbl.(metric);
errorLevelMetric = mean(y_error_metric_values);

fprintf('\nSummary statistics for Error %s:\n', metric);
fprintf('  Mean (used as error level): %.4g\n', errorLevelMetric);
fprintf('  Median: %.4g\n', median(y_error_metric_values, 'omitnan'));
fprintf('  Std:    %.4g\n', std(y_error_metric_values, 'omitnan'));

% --- Get metric display properties ---
[metricLabel, metricUnit] = getMetricDisplayProperties(metric);

% --- Plotting Setup ---
fig = figure('Name', sprintf('%s Trend vs. Error Level', metric), 'Color', 'w', 'Position', [100 100 900 600]);
hold on;

% Plot scatter for trend data
uniqueDepths = unique(x_trend);
colors = lines(numel(uniqueDepths));
% Plot all trend data points with a single legend entry
scatter(x_trend, y_trend_metric, 36, 'filled', 'DisplayName', 'Trend Data');

% --- Exponential Fit for Trend Data: y = a*exp(b*x) ---
fitMask = y_trend_metric > 0 & ~isnan(y_trend_metric) & ~isnan(x_trend);
xFit = x_trend(fitMask);
yFit = y_trend_metric(fitMask);

a = NaN; b = NaN; Rsq = NaN;
x_intersect = NaN; y_intersect = NaN;

if numel(xFit) >= 2
    try
        p = polyfit(xFit, log(yFit), 1); % log(y) = p(1)*x + p(2)
        a = exp(p(2));
        b = p(1);

        % Plot fit over the range of x_trend data
        xFitLine = linspace(min(x_trend), max(x_trend), 200);
        yExpFitLine = a * exp(b * xFitLine);
        plot(xFitLine, yExpFitLine, 'r--', 'LineWidth', 2, 'DisplayName',sprintf('Exp fit (y=%.2ge^{%.2gx})', a, b));

        % R^2 calculation using all points used for fitting
        yExpFitAll = a * exp(b * xFit);
        SSres = sum((yFit - yExpFitAll).^2);
        SStot = sum((yFit - mean(yFit)).^2);
        if SStot > 0 % Avoid division by zero if all yFit values are the same
            Rsq = 1 - SSres/SStot;
        else
            Rsq = 1; % Perfect fit if all points are the same and on the line
        end
        
        fprintf('\nExponential Fit (y = a*exp(b*x)):\n');
        fprintf('  a = %.4g\n', a);
        fprintf('  b = %.4g\n', b);
        fprintf('  R^2 = %.3f\n', Rsq);

    catch ME
        warning(ME.identifier, 'Could not perform exponential fit: %s', ME.message);
        a = NaN; b = NaN; Rsq = NaN; % Ensure these are NaN if fit fails
    end
else
    warning('Not enough data points (found %d, need at least 2) for exponential fit after filtering.', numel(xFit));
end

% --- Plot Error Level ---
x_limits = xlim; % Get current x-axis limits from scatter/fit
if isempty(x_limits) || all(isnan(x_limits)) || x_limits(1) == x_limits(2) % if xlim is not yet set or invalid
    if ~isempty(x_trend)
      x_limits = [min(x_trend)-0.5, max(x_trend)+0.5];
    else % if no trend data at all
      x_limits = [0, 10]; % Default if no trend data
    end
end
plot(x_limits, [errorLevelMetric, errorLevelMetric], 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Error Level (Mean %s = %.2g)', metric, errorLevelMetric));
xlim(x_limits); % Ensure xlims are set for intersection calculation range

% --- Calculate and Annotate Intersection ---
if ~isnan(a) && ~isnan(b) && b ~= 0 % Ensure fit was successful and b is not zero
    % Solve a*exp(b*x) = errorLevelMetric
    % exp(b*x) = errorLevelMetric / a
    % b*x = log(errorLevelMetric / a)
    % x = log(errorLevelMetric / a) / b
    if errorLevelMetric > 0 && a > 0 % Arguments to log must be positive
        x_intersect = log(errorLevelMetric / a) / b;
        y_intersect = errorLevelMetric; % or a * exp(b * x_intersect)

        fprintf('\nIntersection Point:\n');
        fprintf('  MinDiffDepth (x) = %.4g\n', x_intersect);
        fprintf('  %s (y)      = %.4g\n', metric, y_intersect);

        % Plot intersection point if it's within the plotted x-range
        current_xlim = xlim;
        if x_intersect >= current_xlim(1) && x_intersect <= current_xlim(2)
            plot(x_intersect, y_intersect, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Intersection');
            text(x_intersect + 0.05*diff(current_xlim), y_intersect, ...
                 sprintf('Intersection\n(%.2f, %.2g)', x_intersect, y_intersect), ...
                 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'bottom');
        else
            fprintf('  Intersection point (%.4g, %.4g) is outside the current plot x-range [%.2f, %.2f].\n', x_intersect, y_intersect, current_xlim(1), current_xlim(2));
        end
    else
        fprintf('  Cannot calculate intersection: errorLevelMetric / a is not positive (%.2g / %.2g).\n', errorLevelMetric, a);
    end
elseif ~isnan(a) && ~isnan(b) && b == 0 % Case: exponential fit is a horizontal line y = a
     if abs(a - errorLevelMetric) < 1e-9 % Effectively a == errorLevelMetric
        fprintf('\nIntersection: Exponential fit (y=%.2g) is essentially the same as the error level (y=%.2g).\nLines overlap.\n', a, errorLevelMetric);
        % No single intersection point to plot in this case, or infinite.
    else
        fprintf('\nNo intersection: Exponential fit is a horizontal line y=%.2g, and error level is y=%.2g.\nLines are parallel.\n', a, errorLevelMetric);
    end
else
    fprintf('\nCannot calculate intersection because exponential fit was not successful.\n');
end

% --- Finalize Plot ---
hold off;
xlabel('# of Branches Deep');
ylabel(sprintf('%s%s', metricLabel, metricUnit));
title('Effect of Wire on SSTDR vs Branch Depth');
grid on;
legend('show', 'Location', 'best');

% --- Save Figure ---
[~, trendBaseName, ~] = fileparts(trendCsvPath);
[~, errorBaseName, ~] = fileparts(errorCsvPath);
% Use current directory if paths are just filenames
outDirTrend = fileparts(trendCsvPath); 
if isempty(outDirTrend), outDirTrend = pwd; end

outFile = fullfile(outDirTrend, sprintf('%s_vs_%s_%s_intersection.png', trendBaseName, errorBaseName, metric));
try
    saveas(fig, outFile);
    fprintf('Saved plot to %s\n', outFile);
catch ME_save
    fprintf('Could not save plot to %s: %s\nAttempting to save in current directory.\n', outFile, ME_save.message);
    altOutFile = fullfile(pwd, sprintf('%s_vs_%s_%s_intersection.png', trendBaseName, errorBaseName, metric));
    try
        saveas(fig, altOutFile);
        fprintf('Saved plot to %s\n', altOutFile);
    catch ME_save_alt
        fprintf('Failed to save plot in current directory as well: %s\n', ME_save_alt.message);
    end
end

end

% --- Helper Function for Metric Display Properties ---
function [label, unit] = getMetricDisplayProperties(metric)
    switch metric
        case 'PeakMag'
            label = 'Peak Change in Magnitude';
            unit = '';
        case 'AreaNorm'
            label = 'Normalized Area Under Curve';
            unit = '';
        case 'AreaSquared'
            label = 'Squared Area Under Curve';
            unit = '';
        otherwise
            label = metric;
            unit = '';
    end
end 