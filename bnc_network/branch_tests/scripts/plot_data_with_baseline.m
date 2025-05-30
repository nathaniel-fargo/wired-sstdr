% plot_data_with_baseline.m
%
% Plots PeakMag vs. MinDiffDepth from a "trend" CSV, fits an exponential curve,
% plots an error level based on PeakMag from an "error" CSV,
% and calculates and annotates their intersection.
%
% Usage:
%   plot_data_with_baseline('path/to/trend_analysis.csv', 'path/to/error_analysis.csv');
%   plot_data_with_baseline(); % prompts for files

function plot_data_with_baseline(trendCsvPath, errorCsvPath)

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

assert(isfile(trendCsvPath), 'Trend CSV file not found: %s', trendCsvPath);
assert(isfile(errorCsvPath), 'Error CSV file not found: %s', errorCsvPath);

% --- Load and Process Trend Data ---
trendTbl = readtable(trendCsvPath);
fprintf('Loaded trend data from: %s\n', trendCsvPath);

% Remove rows with NaN in relevant columns for trend data
trendRelevantCols = {'MinDiffDepth','PeakMag'};
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
    error('No valid data rows found in trend CSV after NaN removal for MinDiffDepth and PeakMag.');
end

x_trend = trendTbl.MinDiffDepth;
y_trend_peakmag = trendTbl.PeakMag;

fprintf('\nSummary statistics for Trend PeakMag:\n');
fprintf('  Mean:   %.4g\n', mean(y_trend_peakmag, 'omitnan'));
fprintf('  Median: %.4g\n', median(y_trend_peakmag, 'omitnan'));
fprintf('  Std:    %.4g\n', std(y_trend_peakmag, 'omitnan'));

% --- Load and Process Error Data ---
errorTbl = readtable(errorCsvPath);
fprintf('Loaded error data from: %s\n', errorCsvPath);

% Remove rows with NaN in PeakMag for error data
if ~ismember('PeakMag', errorTbl.Properties.VariableNames)
    error('Column "PeakMag" not found in error CSV: %s', errorCsvPath);
end
errorRowMask = ~isnan(errorTbl.PeakMag);
errorTbl = errorTbl(errorRowMask, :);

if isempty(errorTbl)
    error('No valid data rows found in error CSV after NaN removal for PeakMag.');
end

y_error_peakmag_values = errorTbl.PeakMag;
errorLevelPeakMag = mean(y_error_peakmag_values);

fprintf('\nSummary statistics for Error PeakMag:\n');
fprintf('  Mean (used as error level): %.4g\n', errorLevelPeakMag);
fprintf('  Median: %.4g\n', median(y_error_peakmag_values, 'omitnan'));
fprintf('  Std:    %.4g\n', std(y_error_peakmag_values, 'omitnan'));


% --- Plotting Setup ---
fig = figure('Name','PeakMag Trend vs. Error Level','Color','w','Position',[100 100 900 600]);
hold on;

% Plot scatter for trend data
uniqueDepths = unique(x_trend);
colors = lines(numel(uniqueDepths));
% Plot all trend data points with a single legend entry
scatter(x_trend, y_trend_peakmag, 36, 'filled', 'DisplayName', 'Trend Data');
% for idx = 1:numel(uniqueDepths)
%     d = uniqueDepths(idx);
%     scatter(x_trend(x_trend==d), y_trend_peakmag(x_trend==d), 36, colors(idx,:), 'filled', 'DisplayName', sprintf('Trend Data (Depth %d)', d));
% end

% --- Exponential Fit for Trend Data: y = a*exp(b*x) ---
fitMask = y_trend_peakmag > 0 & ~isnan(y_trend_peakmag) & ~isnan(x_trend);
xFit = x_trend(fitMask);
yFit = y_trend_peakmag(fitMask);

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
plot(x_limits, [errorLevelPeakMag, errorLevelPeakMag], 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Error Level (Mean PeakMag = %.2g)', errorLevelPeakMag));
xlim(x_limits); % Ensure xlims are set for intersection calculation range

% --- Calculate and Annotate Intersection ---
if ~isnan(a) && ~isnan(b) && b ~= 0 % Ensure fit was successful and b is not zero
    % Solve a*exp(b*x) = errorLevelPeakMag
    % exp(b*x) = errorLevelPeakMag / a
    % b*x = log(errorLevelPeakMag / a)
    % x = log(errorLevelPeakMag / a) / b
    if errorLevelPeakMag > 0 && a > 0 % Arguments to log must be positive
        x_intersect = log(errorLevelPeakMag / a) / b;
        y_intersect = errorLevelPeakMag; % or a * exp(b * x_intersect)

        fprintf('\nIntersection Point:\n');
        fprintf('  MinDiffDepth (x) = %.4g\n', x_intersect);
        fprintf('  PeakMag (y)      = %.4g\n', y_intersect);

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
        fprintf('  Cannot calculate intersection: errorLevelPeakMag / a is not positive (%.2g / %.2g).\n', errorLevelPeakMag, a);
    end
elseif ~isnan(a) && ~isnan(b) && b == 0 % Case: exponential fit is a horizontal line y = a
     if abs(a - errorLevelPeakMag) < 1e-9 % Effectively a == errorLevelPeakMag
        fprintf('\nIntersection: Exponential fit (y=%.2g) is essentially the same as the error level (y=%.2g).\nLines overlap.\n', a, errorLevelPeakMag);
        % No single intersection point to plot in this case, or infinite.
    else
        fprintf('\nNo intersection: Exponential fit is a horizontal line y=%.2g, and error level is y=%.2g.\nLines are parallel.\n', a, errorLevelPeakMag);
    end
else
    fprintf('\nCannot calculate intersection because exponential fit was not successful.\n');
end

% --- Finalize Plot ---
hold off;
xlabel('# of Branches Deep');
ylabel('Peak Change in Magnitude');
title('Effect of Wire on SSTDR vs Branch Depth');
grid on;
legend('show', 'Location', 'best');

% --- Save Figure ---
[~, trendBaseName, ~] = fileparts(trendCsvPath);
[~, errorBaseName, ~] = fileparts(errorCsvPath);
% Use current directory if paths are just filenames
outDirTrend = fileparts(trendCsvPath); 
if isempty(outDirTrend), outDirTrend = pwd; end

outFile = fullfile(outDirTrend, sprintf('%s_vs_%s_PeakMag_intersection.png', trendBaseName, errorBaseName));
try
    saveas(fig, outFile);
    fprintf('Saved plot to %s\n', outFile);
catch ME_save
    fprintf('Could not save plot to %s: %s\nAttempting to save in current directory.\n', outFile, ME_save.message);
    altOutFile = fullfile(pwd, sprintf('%s_vs_%s_PeakMag_intersection.png', trendBaseName, errorBaseName));
    try
        saveas(fig, altOutFile);
        fprintf('Saved plot to %s\n', altOutFile);
    catch ME_save_alt
        fprintf('Failed to save plot in current directory as well: %s\n', ME_save_alt.message);
    end
end

end 