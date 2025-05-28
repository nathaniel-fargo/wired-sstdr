% draw_network.m
%
% Visualizes branch networks and cable lengths as a directed graph.
% Reads a network CSV and cable logs, and draws each network's hierarchy and cable lengths.
% Can save plots to a folder or plot to a given axes handle.
%
% Author: Nathaniel Fargo
% Date: 2025-05-22
% Org: U of U WIRED
%
% Usage:
%   draw_network();
%   draw_network('my_networks.csv', {'T0','T2'}, 'plots/');
%   draw_network(networks_csv, network_ids, output_folder, ax_handle);
%

function fig_handle = draw_network(networks_csv, network_ids, output_folder, ax_handle)
    % Handle arguments and locate data files ----------------------------------
    if nargin < 1 || isempty(networks_csv)
        thisFile = mfilename('fullpath');
        [thisDir, ~] = fileparts(thisFile);                % .../branch_tests/scripts
        branchTestsDir   = fileparts(thisDir);             % .../branch_tests
        networks_csv = fullfile(branchTestsDir, '2025-05-21', 'networks.csv');
    else
        thisFile = mfilename('fullpath');
        [thisDir, ~] = fileparts(thisFile);                % .../branch_tests/scripts
        branchTestsDir   = fileparts(thisDir);             % .../branch_tests
    end
    if nargin < 2
        network_ids = [];
    end
    if nargin < 3
        output_folder = '';
    end
    if nargin < 4
        ax_handle = []; % Default to no axes handle
    end

    % Cable logs path is always relative to script location
    projectRootDir = fileparts(branchTestsDir);           % .../bnc_network
    logsPath = fullfile(projectRootDir, 'cable_measurements', 'cable_logs.csv');

    assert(isfile(networks_csv), 'Could not find networks.csv at %s', networks_csv);
    assert(isfile(logsPath),     'Could not find cable_logs.csv at %s', logsPath);

    % -------------------------------------------------------------------------
    % Read CSV files ----------------------------------------------------------
    netTbl   = readtable(networks_csv,  'Delimiter', ',');
    logsTbl  = readtable(logsPath,      'Delimiter', ',');

    % Build a map: wire ID -> length string ----------------------------------
    lenMap = containers.Map('KeyType','char','ValueType','char');
    for k = 1:height(logsTbl)
        wire = strtrim(logsTbl.Cable{k});   % e.g. 'A00'
        len  = strtrim(logsTbl.Length{k});  % e.g. 4'
        if ~isempty(wire)
            lenMap(wire) = len;
        end
    end

    % -------------------------------------------------------------------------
    % Filter networks if network_ids provided ---------------------------------
    varNames = netTbl.Properties.VariableNames;
    all_ids = netTbl{:, varNames{1}};
    if isempty(network_ids)
        idxs = 1:height(netTbl);
    else
        idxs = find(ismember(all_ids, network_ids));
    end
    if isempty(idxs)
        warning('No matching network IDs found.');
        fig_handle = []; % Return empty if no networks to plot
        return;
    end

    % Create output folder if needed ------------------------------------------
    if ~isempty(output_folder)
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
    end

    % -------------------------------------------------------------------------
    % Iterate over each selected network and draw -----------------------------
    fig_handle = []; % Initialize fig_handle, will be set to the last figure created or parent of ax_handle

    for ni = 1:numel(idxs)
        n = idxs(ni);
        try
            id  = netTbl{n, varNames{1}}{1};
        catch
            id  = sprintf('#%d', n);
        end
        netStr = netTbl{n, varNames{2}}{1};

        fprintf('Parsing network %s...\n', id);

        [nodeList, edgePairs, terminations] = parse_network_string(netStr);

        % Build parent->children map for traversal
        parentMap = containers.Map('KeyType','char','ValueType','any');
        childSet = containers.Map('KeyType','char','ValueType','logical');
        for i = 1:size(edgePairs,1)
            p = edgePairs{i,1};
            c = edgePairs{i,2};
            if ~isKey(parentMap, p)
                parentMap(p) = {};
            end
            children = parentMap(p);
            children{end+1} = c;
            parentMap(p) = children;
            childSet(c) = true;
        end
        % Find root (wire that is never a child)
        root = '';
        for i = 1:numel(nodeList)
            w = nodeList{i};
            if ~isKey(childSet, w)
                root = w;
                break;
            end
        end
        if isempty(root)
            root = nodeList{1}; % fallback
        end

        % Plot setup
        current_ax = [];
        if ~isempty(ax_handle) && isgraphics(ax_handle, 'axes')
            axes(ax_handle); % Make the provided axes current
            current_ax = ax_handle;
            fig_handle = get(ax_handle, 'Parent');
            % For plotting on existing axes, usually only one network is drawn.
            % If multiple are in idxs, this will overwrite. Caller should manage this.
            if numel(idxs) > 1
                warning('Multiple network IDs provided with a single axes handle. Only the last network will be visible on the axes.');
            end
        else
            hFig = figure('Name', ['Network ' char(id)], 'NumberTitle', 'off');
            fig_handle = hFig; % Store handle of the created figure
            current_ax = gca; % Get current axes of the new figure
        end
        
        hold(current_ax, 'on');
        axis(current_ax, 'equal', 'off');
        title(current_ax, ['Network ' char(id)]);

        % Recursively draw the network as lines (left to right)
        baseAngle = 0; % horizontal right
        baseLen = 1.5;    % base length for each segment
        draw_branch_recursive(current_ax, root, [0,0], baseAngle, 0, parentMap, terminations, lenMap, 0, baseLen);
        hold(current_ax, 'off');

        % Save if requested and if we created the figure
        if ~isempty(output_folder) && isempty(ax_handle)
            saveas(fig_handle, fullfile(output_folder, sprintf('network_%s.png', char(id))));
            % Do not close here if we are in a loop and returning the last fig_handle
            % The original script also did close all at the start.
        end
    end

    % If not plotting to a specific ax_handle and no output folder, figures remain open.
    % If plotting to a specific ax_handle, caller manages figure visibility.
    % If output_folder is specified, figures are saved. Original script had close all at the start.

    end

    function draw_text_with_bg(ax, x, y, str, varargin)
    % Draw text with a fixed-size white rectangle background for readability
        rect_w = 2; % fixed width
        rect_h = 1; % fixed height
        % Draw white rectangle
        rectangle(ax, 'Position', [x-rect_w/2, y-rect_h/2, rect_w, rect_h], ...
            'FaceColor', 'w', 'EdgeColor', 'none', 'Curvature', 0.1);
        % Draw text centered
        text(ax, x, y, str, varargin{:}, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end

    function draw_branch_recursive(ax, wire, startPt, angle, depth, parentMap, termMap, lenMap, branchIdx, baseLen)
    % Recursively draw a wire as a line, with junctions and terminations
        % Determine length (for plotting, not to scale)
        if lenMap.isKey(wire)
            lenStr = lenMap(wire);
            % Try to extract feet as a number
            ft = regexp(lenStr, '(\d+)', 'match');
            if ~isempty(ft)
                segLen = str2double(ft{1}) * 0.2 + 1.2; % scale for visual
            else
                segLen = baseLen;
            end
        else
            lenStr = '??';
            segLen = baseLen;
        end
        % Compute end point
        dx = segLen * cos(angle);
        dy = segLen * sin(angle);
        endPt = startPt + [dx, dy];

        % Draw the wire as a line
        plot(ax, [startPt(1), endPt(1)], [startPt(2), endPt(2)], 'k-', 'LineWidth', 2);
        % Draw label at midpoint, offset vertically for left-right layout
        midPt = (startPt + endPt)/2;
        labelOffset = [0, 0.25]; % more vertical space
        draw_text_with_bg(ax, midPt(1)+labelOffset(1), midPt(2)+labelOffset(2), sprintf('%s\n%s', wire, lenStr), ...
            'FontSize', 9);

        % If this wire has children, draw a junction dot at end
        if isKey(parentMap, wire)
            plot(ax, endPt(1), endPt(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
            nChild = numel(parentMap(wire));
            if nChild == 1
                childAngles = angle;
            elseif nChild == 2
                childAngles = [-pi/2, 0]; % first child down, second right
            else
                spread = pi/2; % fallback: spread angle for branches (vertical)
                childAngles = linspace(angle-spread/2, angle+spread/2, nChild);
            end
            children = parentMap(wire);
            for c = 1:nChild
                draw_branch_recursive(ax, children{c}, endPt, childAngles(c), depth+1, parentMap, termMap, lenMap, c, baseLen*0.95);
            end
        else
            % Termination: draw a fixed-size red rectangle with white text at end, no brackets
            if termMap.isKey(wire)
                t = termMap(wire); % just the letter, no brackets
            else
                t = '';
            end
            termOffset = [0.25, 0]; % more horizontal space
            rect_w = 0.7; % fixed width for termination
            rect_h = 0.35; % fixed height for termination
            % Draw red rectangle
            rectangle(ax, 'Position', [endPt(1)+termOffset(1)-rect_w/2, endPt(2)+termOffset(2)-rect_h/2, rect_w, rect_h], ...
                'FaceColor', 'r', 'EdgeColor', 'none', 'Curvature', 0.1);
            % Draw the white text on top, centered
            text(ax, endPt(1)+termOffset(1), endPt(2)+termOffset(2), t, ...
                'FontSize', 10, 'FontWeight','bold', 'Color','w', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
        end
    end

    % =========================================================================
    function [nodes, edges, termMap] = parse_network_string(str)
    % PARSE_NETWORK_STRING  Convert network definition into node/edge lists.
    %   [nodes, edges, termMap] = parse_network_string(str)
    %   * nodes   : 1×N cell array of unique wire IDs
    %   * edges   : M×2 cell array, each row = {parent, child}
    %   * termMap : containers.Map with termination code for leaf wires

    str = regexprep(str, '\s', '');       % remove all whitespace

    nodes = {};
    edges = {};
    termMap = containers.Map('KeyType','char','ValueType','char');
    pos = 1;
    len = numel(str);

        function [wire] = parseNode()
            %% Expect opening '{'
            if str(pos) ~= '{'
                error('Expected "{" at position %d in: %s', pos, str);
            end
            pos = pos + 1;   % consume '{'

            % ---- Read wire name --------------------------------------------
            start = pos;
            while pos <= len && ~ismember(str(pos), ['[', '{', '}', ' '])
                pos = pos + 1;
            end
            wire = str(start:pos-1);
            if isempty(wire)
                error('Empty wire name at position %d', pos);
            end
            if ~any(strcmp(nodes, wire))
                nodes{end+1} = wire; %#ok<AGROW>
            end

            % ---- Optional termination --------------------------------------
            if pos <= len && str(pos) == '['
                pos = pos + 1; % consume '['
                tStart = pos;
                while pos <= len && str(pos) ~= ']'
                    pos = pos + 1;
                end
                if pos > len
                    error('Unclosed termination bracket for wire %s', wire);
                end
                termCode = str(tStart:pos-1);
                termMap(wire) = termCode;
                pos = pos + 1; % consume ']'
            end

            % ---- Parse zero or more child branches -------------------------
            while pos <= len && str(pos) == '{'
                child = parseNode();
                edges(end+1, :) = {wire, child}; %#ok<AGROW>
            end

            % ---- Expect closing '}' ----------------------------------------
            if pos > len || str(pos) ~= '}'
                error('Expected "}" for wire %s (position %d)', wire, pos);
            end
            pos = pos + 1; % consume '}'
        end

    % Kick-off parsing --------------------------------------------------------
    parseNode();

    % Verify full consumption -------------------------------------------------
    if pos <= len
        remaining = str(pos:end);
        warning('Parser did not consume entire string. Remaining: %s', remaining);
    end
    end
