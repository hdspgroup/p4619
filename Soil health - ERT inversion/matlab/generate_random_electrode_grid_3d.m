clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));

outDir = fullfile(rootDir, 'data', 'electrode_arrays');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Configuration

% rngSeed = 20260512;
rngSeed = 20260521;
nElectrodes = 45;
maxRadius = 25;
centerXY = [0, 0];
electrodeZ = 0;
currentA = 1;
targetDepth = 3;
progressMode = "text"; % text | waitbar | off

% Practical quadrupole filters for a random 2D surface array.
options = struct();
options.MinCurrentSpacing = 2.5;
options.MinPotentialSpacing = 2.5;
options.MinCentroidSeparation = 2.5;
options.MaxCurrentSpacing = 50;
options.MaxPotentialSpacing = 50;
options.MaxGeometricFactor = 2e4;
options.MinSensitivity = 1e-4;
options.TargetDepth = targetDepth;
options.DepthToCurrentSpacingFactor = 3.0;
options.DepthToPotentialSpacingFactor = 2.0;
options.DepthToCentroidSpacingFactor = 3.0;

arrayTag = sprintf('random_planar_n%d_r%d_seed%d', nElectrodes, round(maxRadius), rngSeed);

%% Electrode coordinates

electrodes = generate_random_planar_electrodes( ...
    'NumElectrodes', nElectrodes, ...
    'MaxRadius', maxRadius, ...
    'CenterXY', centerXY, ...
    'Z', electrodeZ, ...
    'Seed', rngSeed);

%% Survey / ABMN combinations

survey = make_ert_planar_survey( ...
    'Electrodes', electrodes, ...
    'Current', currentA, ...
    'MinCurrentSpacing', options.MinCurrentSpacing, ...
    'MinPotentialSpacing', options.MinPotentialSpacing, ...
    'MinCentroidSeparation', options.MinCentroidSeparation, ...
    'MaxCurrentSpacing', options.MaxCurrentSpacing, ...
    'MaxPotentialSpacing', options.MaxPotentialSpacing, ...
    'MaxGeometricFactor', options.MaxGeometricFactor, ...
    'MinSensitivity', options.MinSensitivity, ...
    'ProgressMode', progressMode, ...
    'TargetDepth', options.TargetDepth, ...
    'DepthToCurrentSpacingFactor', options.DepthToCurrentSpacingFactor, ...
    'DepthToPotentialSpacingFactor', options.DepthToPotentialSpacingFactor, ...
    'DepthToCentroidSpacingFactor', options.DepthToCentroidSpacingFactor);

summary = struct();
summary.arrayTag = arrayTag;
summary.rngSeed = rngSeed;
summary.nElectrodes = survey.nElectrodes;
summary.nMeasurements = survey.nMeasurements;
summary.maxRadius = maxRadius;
summary.currentA = currentA;
summary.targetDepth = targetDepth;
summary.progressMode = progressMode;
summary.options = options;

%% Plot

fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08, 0.12, 0.80, 0.66]);
t = linspace(0, 2*pi, 300);

ax1 = subplot(1, 2, 1);
hold(ax1, 'on');
plot(ax1, centerXY(1) + maxRadius * cos(t), centerXY(2) + maxRadius * sin(t), 'k--', 'LineWidth', 1.2);
scatter(ax1, electrodes(:, 1), electrodes(:, 2), 55, 'filled', 'MarkerFaceColor', [0.1 0.35 0.85]);
text(ax1, electrodes(:, 1) + 0.5, electrodes(:, 2) + 0.5, string(1:nElectrodes), 'FontSize', 8, 'Color', [0.15 0.15 0.15]);
plot(ax1, centerXY(1), centerXY(2), 'k+', 'MarkerSize', 12, 'LineWidth', 1.2);
axis(ax1, 'equal');
grid(ax1, 'on');
xlabel(ax1, 'x (m)');
ylabel(ax1, 'y (m)');
title(ax1, sprintf('Electrode layout: %d electrodes', survey.nElectrodes));

ax2 = subplot(1, 2, 2);
hold(ax2, 'on');
plot(ax2, centerXY(1) + maxRadius * cos(t), centerXY(2) + maxRadius * sin(t), 'k--', 'LineWidth', 1.2);
localPlotPairGraph(ax2, electrodes, survey.pairUsage.potential, [0.72 0.82 0.94], 0.5, 2.0);
localPlotPairGraph(ax2, electrodes, survey.pairUsage.current, [0.18 0.42 0.76], 0.8, 2.6);
usage = survey.electrodeUsage.total(:);
markerSizes = 35 + 90 * sqrt(usage ./ max(max(usage), 1));
scatter(ax2, electrodes(:, 1), electrodes(:, 2), markerSizes, usage, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
unusedMask = usage == 0;
if any(unusedMask)
    scatter(ax2, electrodes(unusedMask, 1), electrodes(unusedMask, 2), 120, 'o', ...
        'MarkerEdgeColor', [0.85 0.1 0.1], 'LineWidth', 1.4);
end
text(ax2, electrodes(:, 1) + 0.5, electrodes(:, 2) + 0.5, string(1:nElectrodes), 'FontSize', 8, 'Color', [0.15 0.15 0.15]);
plot(ax2, centerXY(1), centerXY(2), 'k+', 'MarkerSize', 12, 'LineWidth', 1.2);
axis(ax2, 'equal');
grid(ax2, 'on');
xlabel(ax2, 'x (m)');
ylabel(ax2, 'y (m)');
title(ax2, sprintf('Coverage: %d ABMN, %d unused electrodes', ...
    survey.nMeasurements, survey.summary.unusedElectrodes));
cb = colorbar(ax2);
cb.Label.String = 'electrode usage count';
colormap(ax2, parula);
legend(ax2, {'radius bound', 'potential pairs', 'current pairs', 'used electrodes', 'unused electrodes'}, ...
    'Location', 'southoutside');

%% Save

save(fullfile(outDir, [arrayTag '.mat']), 'survey', 'summary');
exportgraphics(fig, fullfile(outDir, [arrayTag '.png']), 'Resolution', 220);

electrodeTable = array2table(electrodes, 'VariableNames', {'x', 'y', 'z'});
writetable(electrodeTable, fullfile(outDir, [arrayTag '_electrodes.csv']));

quadTable = array2table(survey.quads, 'VariableNames', {'A', 'B', 'M', 'N'});
quadTable.K = survey.geometryFactorK(:);
quadTable.geomTerm = survey.geometryTerm(:);
quadTable.distAB = survey.distanceAB(:);
quadTable.distMN = survey.distanceMN(:);
quadTable.centroidSep = survey.centroidSeparation(:);
writetable(quadTable, fullfile(outDir, [arrayTag '_quads.csv']));

function localPlotPairGraph(ax, electrodes, pairMatrix, colorBase, minWidth, maxWidth)
    pairMatrix = triu(pairMatrix, 1);
    [ii, jj, vv] = find(pairMatrix);
    if isempty(vv)
        return;
    end

    vScaled = vv ./ max(vv);
    for k = 1:numel(vv)
        p1 = electrodes(ii(k), 1:2);
        p2 = electrodes(jj(k), 1:2);
        lineColor = 1 - (1 - colorBase) * (0.35 + 0.65 * vScaled(k));
        line(ax, [p1(1), p2(1)], [p1(2), p2(2)], ...
            'Color', lineColor, ...
            'LineWidth', minWidth + (maxWidth - minWidth) * vScaled(k));
    end
end
