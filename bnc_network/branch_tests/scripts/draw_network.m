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
%   draw_network(networks_csv, network_ids, output_folder, ax_handle, wire_labels);
%

function fig_handle = draw_network(networks_csv, network_ids, output_folder, ax_handle, wire_labels)
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
    if nargin < 5 || isempty(wire_labels)
        wire_labels = true; % Default to showing wire labels
    end

    % Cable logs path is always relative to script location
    projectRootDir = fileparts(branchTestsDir);           % .../bnc_network
    logsPath = fullfile(projectRootDir, 'cable_measurements', 'cable_logs.csv');

    assert(isfile(networks_csv), 'Could not find networks.csv at %s', networks_csv);
    assert(isfile(logsPath),     'Could not find cable_logs.csv at %s', logsPath);

    % -------------------------------------------------------------------------
    % Read CSV files ----------------------------------------------------------
    netTbl   = readtable(networks_csv,  'Delimiter', ';', 'TextType', 'string');
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
            % Extract ID - handle both cell arrays and string arrays
            if iscell(netTbl{n, varNames{1}})
                id = netTbl{n, varNames{1}}{1};
            else
                id = netTbl{n, varNames{1}};
            end
        catch
            id  = sprintf('#%d', n);
        end
        
        % Extract network string - handle both cell arrays and string arrays
        if iscell(netTbl{n, varNames{2}})
            netStr = netTbl{n, varNames{2}}{1};
        else
            netStr = netTbl{n, varNames{2}};
        end
        
        % Convert to char array for consistency
        netStr = char(netStr);
        netStr = strtrim(netStr);  % Remove leading/trailing whitespace

        fprintf('Parsing network %s: %s...\n', char(id), netStr);

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
        draw_branch_recursive(current_ax, root, [0,0], baseAngle, 0, parentMap, terminations, lenMap, 0, baseLen, wire_labels);
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
    % Draw text (background removed)
        % Rectangle drawing removed
        % Draw text centered
        text(ax, x, y, str, varargin{:}, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end

    function draw_branch_recursive(ax, wire, startPt, angle, depth, parentMap, termMap, lenMap, branchIdx, baseLen, show_labels)
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

        % Draw the wire as a line (now nice blue, thinner)
        plot(ax, [startPt(1), endPt(1)], [startPt(2), endPt(2)], '-', 'Color', [0.2 0.5 0.8], 'LineWidth', 1.5);
        
        % Draw label at midpoint only if show_labels is true
        if show_labels
            midPt = (startPt + endPt)/2;
            labelOffset = [0, 0.25]; % more vertical space
            draw_text_with_bg(ax, midPt(1)+labelOffset(1), midPt(2)+labelOffset(2), sprintf('%s\n%s', wire, lenStr), ...
                'FontSize', 9);
        end

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
                draw_branch_recursive(ax, children{c}, endPt, childAngles(c), depth+1, parentMap, termMap, lenMap, c, baseLen*0.95, show_labels);
            end
        else
            % Termination: Draw a shape or text based on type
            term_text = '';
            if termMap.isKey(wire)
                term_text = strtrim(termMap(wire));
            end

            markerSize = 8;
            shapeColor = 'r'; % Red for all termination markers/text

            switch lower(term_text)
                case 'open'
                    plot(ax, endPt(1), endPt(2), 'o', 'MarkerEdgeColor', shapeColor, 'MarkerSize', markerSize, 'LineWidth', 1.5);
                case 'short'
                    plot(ax, endPt(1), endPt(2), 'o', 'MarkerEdgeColor', shapeColor, 'MarkerFaceColor', shapeColor, 'MarkerSize', markerSize);
                case 'terminated'
                    plot(ax, endPt(1), endPt(2), 's', 'MarkerEdgeColor', shapeColor, 'MarkerFaceColor', shapeColor, 'MarkerSize', markerSize);
                otherwise
                    % Draw the termination text in red font, centered
                    text(ax, endPt(1), endPt(2), term_text, ...
                        'FontSize', 10, 'FontWeight','bold', 'Color',shapeColor, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
            end
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
            while pos <= len && ~ismember(str(pos), ['[', '{', '}', ',', ' '])
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

            % ---- Parse child branches (comma-separated within braces) -----
            if pos <= len && str(pos) == '{'
                pos = pos + 1; % consume opening '{'
                
                % Parse comma-separated children until closing '}'
                while pos <= len && str(pos) ~= '}'
                    % Each child in the comma-separated list could be:
                    % 1. Simple wire with termination: WireID[term]
                    % 2. Wire with its own nested children: WireID{...}
                    
                    child = parseCommaListItem();
                    edges(end+1, :) = {wire, child}; %#ok<AGROW>
                    
                    % Check for comma or end
                    if pos <= len && str(pos) == ','
                        pos = pos + 1; % consume ','
                    elseif pos <= len && str(pos) == '}'
                        break; % end of children list
                    else
                        error('Expected "," or "}" after child for wire %s at position %d', wire, pos);
                    end
                end
                
                if pos > len || str(pos) ~= '}'
                    error('Expected closing "}" for children of wire %s', wire);
                end
                pos = pos + 1; % consume closing '}'
            end

            % ---- Expect closing '}' ----------------------------------------
            if pos > len || str(pos) ~= '}'
                error('Expected "}" for wire %s (position %d)', wire, pos);
            end
            pos = pos + 1; % consume '}'
        end
        
        function [wire] = parseCommaListItem()
            % Parse an item in a comma-separated list - could be simple or have nested children
            
            % ---- Read wire name --------------------------------------------
            start = pos;
            while pos <= len && ~ismember(str(pos), ['[', '{', '}', ',', ' '])
                pos = pos + 1;
            end
            wire = str(start:pos-1);
            if isempty(wire)
                error('Empty wire name in comma list at position %d', pos);
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

            % ---- Check for nested children --------------------------------
            if pos <= len && str(pos) == '{'
                % This wire has its own children - parse them in the same way as main parsing
                pos = pos + 1; % consume opening '{'
                
                % Parse comma-separated children until closing '}'
                while pos <= len && str(pos) ~= '}'
                    childWire = parseCommaListItem(); % recursive call for each child
                    edges(end+1, :) = {wire, childWire}; %#ok<AGROW>
                    
                    % Check for comma or end
                    if pos <= len && str(pos) == ','
                        pos = pos + 1; % consume ','
                    elseif pos <= len && str(pos) == '}'
                        break; % end of children list
                    else
                        error('Expected "," or "}" after child for wire %s at position %d', wire, pos);
                    end
                end
                
                if pos > len || str(pos) ~= '}'
                    error('Expected closing "}" for children of wire %s', wire);
                end
                pos = pos + 1; % consume closing '}'
            end
            % If no '{', this is just a simple wire (possibly terminated)
        end

    % Kick-off parsing --------------------------------------------------------
    parseNode();

    % Verify full consumption -------------------------------------------------
    if pos <= len
        remaining = str(pos:end);
        warning('Parser did not consume entire string. Remaining: %s', remaining);
    end
    end
