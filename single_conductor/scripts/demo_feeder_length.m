% demo_feeder_length.m
%
% Demo script showing how to use the new feeder length parameter
% in both plot_livewire_folder_time and plot_folder_difference_time
%
% Author: Nathaniel Fargo
% Date: 2025-01-23

clear; close all; clc;

% Get the directory of this script and add it to path
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

fprintf('=== Feeder Length Demo ===\n\n');

% Demo parameters
folderPath = fullfile('..', '06-16-2025', 'LiveWire', 'CSV');
frequency = '48 MHz';
confidence_level = 0.95;
feeder_length = 30; % feet

% Example 1: Groups plot with feeder length
fprintf('1. Plotting groups with feeder length (l=%g ft)...\n', feeder_length);
fig1 = plot_livewire_folder_time(folderPath, frequency, confidence_level, feeder_length);
saveas(fig1, 'demo_groups_with_feeder.png');
close(fig1);

% Example 2: Groups plot without feeder length
fprintf('2. Plotting groups without feeder length...\n');
fig2 = plot_livewire_folder_time(folderPath, frequency, confidence_level);
saveas(fig2, 'demo_groups_no_feeder.png');
close(fig2);

% Example 3: Differences plot with feeder length
fprintf('3. Plotting differences with feeder length (l=%g ft)...\n', feeder_length);
fig3 = plot_folder_difference_time(folderPath, frequency, confidence_level, feeder_length);
saveas(fig3, 'demo_differences_with_feeder.png');
close(fig3);

% Example 4: Differences plot without feeder length
fprintf('4. Plotting differences without feeder length...\n');
fig4 = plot_folder_difference_time(folderPath, frequency, confidence_level);
saveas(fig4, 'demo_differences_no_feeder.png');
close(fig4);

% Example 5: Multiple feeder lengths comparison
fprintf('5. Comparing different feeder lengths...\n');
feeder_lengths = [20, 30, 50]; % feet

figure;
for i = 1:length(feeder_lengths)
    subplot(1, 3, i);
    
    % Create a temporary plot for each feeder length
    temp_fig = plot_livewire_folder_time(folderPath, frequency, confidence_level, feeder_lengths(i));
    title(sprintf('Feeder Length: %g ft', feeder_lengths(i)));
    
    % Copy the plot to the subplot
    temp_ax = gca(temp_fig);
    copyobj(allchild(temp_ax), gca);
    xlim(temp_ax.XLim);
    ylim(temp_ax.YLim);
    xlabel('Time (ns)');
    ylabel('Normalized Value');
    legend('show', 'Location', 'best');
    grid on;
    
    close(temp_fig);
end

sgtitle('Comparison of Different Feeder Lengths');
saveas(gcf, 'demo_feeder_comparison.png');
close(gcf);

fprintf('\n=== Summary ===\n');
fprintf('New features added:\n');
fprintf('✓ Feeder cable length parameter (l) in feet\n');
fprintf('✓ Vertical dashed line showing feeder length time\n');
fprintf('✓ Updated title format: "...at 48 MHz (l=30ft)"\n');
fprintf('✓ Legend entry: "Feeder Line Length"\n');
fprintf('✓ Control group (Z00) labeled as "Control"\n');
fprintf('✓ Neutral gray color for control groups\n');
fprintf('✓ Backward compatibility when l parameter not provided\n');

% Calculate and display time conversion info
feeder_time_ns = (feeder_length * 2) / (983571056 * 0.66) * 1e9;
fprintf('\nTime conversion example:\n');
fprintf('Feeder length: %g ft\n', feeder_length);
fprintf('Round-trip time: %.2f ns\n', feeder_time_ns);
fprintf('(Using VOP=0.66, speed of light, factor of 2 for round trip)\n');

fprintf('\nDemo completed! Check the generated PNG files.\n'); 