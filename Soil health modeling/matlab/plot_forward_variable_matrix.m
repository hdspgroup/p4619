function plot_forward_variable_matrix()
% Plot a scatter matrix for forward-model variables with color coding:
% - blue when both variables are state variables or both are context variables
% - red when one variable is a state variable and the other is a context variable
%
% Uses the stable linked point table when available so EC and CEC can appear
% in the same table.

    close all force;
    clc;

    paths = codex_paths('forward_variable_matrix');
    stablePath = fullfile(paths.outputsDir, 'ec_cec_linkage', 'stable_point_table.csv');

    if isfile(stablePath)
        T = readtable(stablePath);
    else
        s = load(fullfile(paths.dataDir, 'data_comp.mat'));
        varNames = matlab.lang.makeValidName(string(s.varNames), 'ReplacementStyle', 'delete');
        T = array2table(s.stableData, 'VariableNames', cellstr(varNames));
    end

    if ismember("EC", string(T.Properties.VariableNames))
        T.EC = T.EC .* 0.1; % dS/m to S/m
    end

    stateVars = ["EC","N","P","K","CEC","OC","Moisture"];
    contextVars = ["Clay","Silt","Sand","Coarse","pH_CaCl2","pH_H2O","CaCO3"];
    moistureAliases = ["Moisture","H_2O_Vol_","H_2OVol","H_2O_Vol","WaterContent","theta_v"];

    for a = moistureAliases
        if ismember(a, string(T.Properties.VariableNames)) && a ~= "Moisture"
            T.Moisture = T.(a);
        end
    end

    presentState = intersect(stateVars, string(T.Properties.VariableNames), 'stable');
    presentContext = intersect(contextVars, string(T.Properties.VariableNames), 'stable');
    vars = [presentState, presentContext];

    D = T(:, vars);
    D = D(all(~ismissing(D), 2), :);
    D = D(all(isfinite(D{:,:}), 2), :);

    maxPoints = 2500;
    rng(7);
    if height(D) > maxPoints
        D = D(randperm(height(D), maxPoints), :);
    end

    n = numel(vars);
    stateSet = presentState;
    contextSet = presentContext;
    sameColor = [0.12 0.39 0.82];
    crossColor = [0.82 0.22 0.18];
    histColor = [0.55 0.55 0.60];

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 170*n 170*n]);
    tl = tiledlayout(n, n, 'TileSpacing', 'compact', 'Padding', 'compact');

    for r = 1:n
        for c = 1:n
            ax = nexttile((r - 1) * n + c);
            x = D.(vars(c));
            y = D.(vars(r));

            if r == c
                histogram(ax, x, 20, 'FaceColor', histColor, 'EdgeColor', 'none');
            else
                if (any(vars(r) == stateSet) && any(vars(c) == stateSet)) || ...
                   (any(vars(r) == contextSet) && any(vars(c) == contextSet))
                    thisColor = sameColor;
                else
                    thisColor = crossColor;
                end
                scatter(ax, x, y, 8, 'filled', ...
                    'MarkerFaceColor', thisColor, 'MarkerFaceAlpha', 0.25, ...
                    'MarkerEdgeAlpha', 0.25);
            end

            grid(ax, 'on');
            ax.FontSize = 8;

            if r == 1
                title(ax, localPrettyName(vars(c)), 'Interpreter', 'none', 'FontSize', 9);
            else
                ax.Title.String = '';
            end

            if c == 1
                ylabel(ax, localPrettyName(vars(r)), 'Interpreter', 'none', 'FontSize', 9);
            else
                ax.YTickLabel = [];
            end

            if r < n
                ax.XTickLabel = [];
            else
                xlabel(ax, localPrettyName(vars(c)), 'Interpreter', 'none', 'FontSize', 8);
            end
        end
    end

    % Create invisible objects for a figure-level legend.
    hold(nexttile(1), 'on');
    h1 = scatter(nan, nan, 30, 'filled', 'MarkerFaceColor', sameColor);
    h2 = scatter(nan, nan, 30, 'filled', 'MarkerFaceColor', crossColor);
    lg = legend([h1 h2], {'State-State / Context-Context', 'State-Context'}, ...
        'Orientation', 'horizontal', 'Location', 'southoutside');
    lg.Box = 'off';

    title(tl, 'Forward Model Variable Matrix', 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(f, fullfile(paths.analysisDir, 'forward_variable_matrix.png'), 'Resolution', 220);
    close(f);

    fid = fopen(fullfile(paths.analysisDir, 'forward_variable_matrix_legend.txt'), 'w');
    fprintf(fid, 'Forward variable matrix\n');
    fprintf(fid, 'Blue panels: both variables are in the same block (state-state or context-context).\n');
    fprintf(fid, 'Red panels: one variable is a state variable and the other is a context variable.\n');
    fprintf(fid, 'Variables included: %s\n', strjoin(cellstr(vars), ', '));
    fclose(fid);
end

function label = localPrettyName(name)
    switch char(name)
        case 'pH_CaCl2'
            label = 'pH CaCl2';
        case 'pH_H2O'
            label = 'pH H2O';
        otherwise
            label = char(name);
    end
end
