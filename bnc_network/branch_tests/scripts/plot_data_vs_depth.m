% plot_data_vs_depth.m
%
% Plots PeakMag, AreaUnsquared, and AreaSquared against MinDiffDepth from a network analysis CSV.
% Usage:
%   plot_data_vs_depth('path/to/networks_n_analysis.csv');
%   plot_data_vs_depth(); % prompts for file

function plot_data_vs_depth(csvPath, statsToPlot)

if nargin < 1 || isempty(csvPath)
    [f,p] = uigetfile('*.csv', 'Select analysis CSV');
    if isequal(f,0), disp('No file selected.'); return; end
    csvPath = fullfile(p,f);
end

if nargin < 2 || isempty(statsToPlot)
    statsToPlot = {'PeakMag','AreaUnsquared','AreaSquared'};
end

assert(isfile(csvPath), 'CSV file not found: %s', csvPath);
tbl = readtable(csvPath);

% --- Remove rows with NaN in any relevant column ---
relevantCols = {'MinDiffDepth','PeakMag','AreaUnsquared','AreaSquared'};
rowMask = true(height(tbl),1);
for k = 1:numel(relevantCols)
    if ismember(relevantCols{k}, tbl.Properties.VariableNames)
        rowMask = rowMask & ~isnan(tbl.(relevantCols{k}));
    end
end
tbl = tbl(rowMask, :);

% --- Print summary statistics for selected stats ---
statPrintMap = struct('PeakMag', 'PeakMag', ...
                     'AreaUnsquared', 'AreaUnsquared', ...
                     'AreaSquared', 'AreaSquared');
for s = 1:numel(statsToPlot)
    stat = statsToPlot{s};
    if ~isfield(statPrintMap, stat)
        warning('Unknown statistic: %s', stat);
        continue;
    end
    y = tbl.(statPrintMap.(stat));
    if strcmp(stat, 'AreaUnsquared')
        y = abs(y); % use absolute value for AreaUnsquared
    end
    fprintf('\nSummary statistics for %s:\n', stat);
    fprintf('  Mean:   %.4g\n', mean(y, 'omitnan'));
    fprintf('  Median: %.4g\n', median(y, 'omitnan'));
    fprintf('  Std:    %.4g\n', std(y, 'omitnan'));
end

% Check required columns
reqCols = {'MinDiffDepth','PeakMag','AreaUnsquared','AreaSquared'};
for k = 1:numel(reqCols)
    if ~ismember(reqCols{k}, tbl.Properties.VariableNames)
        error('Column "%s" not found in %s', reqCols{k}, csvPath);
    end
end

% Extract data
x = tbl.MinDiffDepth;
peakMag = tbl.PeakMag;
areaUnsq = abs(tbl.AreaUnsquared); % always use absolute value
areaSq = tbl.AreaSquared;

% Map stat names to data and labels
statMap = struct('PeakMag', peakMag, ...
                 'AreaUnsquared', areaUnsq, ...
                 'AreaSquared', areaSq);
ylabels = struct('PeakMag','PeakMag', ...
                 'AreaUnsquared','|AreaUnsquared|', ...
                 'AreaSquared','AreaSquared');
titles = struct('PeakMag','PeakMag vs. MinDiffDepth', ...
                'AreaUnsquared','|AreaUnsquared| vs. MinDiffDepth', ...
                'AreaSquared','AreaSquared vs. MinDiffDepth');

nStats = numel(statsToPlot);
fig = figure('Name','Statistics vs. MinDiffDepth','Color','w','Position',[100 100 900 300*nStats]);
uniqueDepths = unique(x);
colors = lines(numel(uniqueDepths));

for s = 1:nStats
    stat = statsToPlot{s};
    y = statMap.(stat);
    subplot(nStats,1,s);
    hold on;
    for idx = 1:numel(uniqueDepths)
        d = uniqueDepths(idx);
        scatter(x(x==d), y(x==d), 36, colors(idx,:), 'filled', 'DisplayName', sprintf('Depth %d', d));
    end
    % Exponential fit: y = a*exp(b*x)
    % Use all individual data points where y > 0
    fitMask = y > 0 & ~isnan(y) & ~isnan(x);
    xFit = x(fitMask);
    yFit = y(fitMask);
    if numel(xFit) >= 2
        p = polyfit(xFit, log(yFit), 1); % log(y) = p(1)*x + p(2)
        a = exp(p(2));
        b = p(1);
        % Plot fit over the range of x
        xFitLine = linspace(min(xFit), max(xFit), 100);
        yExpFitLine = a * exp(b * xFitLine);
        plot(xFitLine, yExpFitLine, 'r--', 'LineWidth', 2, 'DisplayName','Exp fit');
        % R^2 calculation using all points
        yExpFitAll = a * exp(b * xFit);
        SSres = sum((yFit - yExpFitAll).^2);
        SStot = sum((yFit - mean(yFit)).^2);
        Rsq = 1 - SSres/SStot;
        % Annotate equation and R^2
        eqnStr = sprintf('y = %.2g e^{%.2g x}\nR^2 = %.3f', a, b, Rsq);
        xText = min(xFit) + 0.05*(max(xFit)-min(xFit));
        yText = max(yFit) - 0.1*(max(yFit)-min(yFit));
        text(xText, yText, eqnStr, 'FontSize', 10, 'Color', 'r', 'VerticalAlignment','top');
    end
    hold off;
    xlabel('MinDiffDepth'); ylabel(ylabels.(stat)); title(titles.(stat)); grid on;
    legend('show');
end

sgtitle('Network Difference Statistics vs. MinDiffDepth');

% Save figure
[outDir, baseName, ~] = fileparts(csvPath);
outFile = fullfile(outDir, baseName + "_stats_vs_depth.png");
saveas(fig, outFile);
fprintf('Saved plot to %s\n', outFile);

end 