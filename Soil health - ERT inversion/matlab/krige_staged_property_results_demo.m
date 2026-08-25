clear; clc; close all;

% Krige staged soil-property inversion results onto high-resolution slices.
% The default is slice kriging, not a full 3-D volume, because local ordinary
% kriging requires solving a small system at every output point.

saveArg = true;
fieldSource = "current"; % current | actual | background | start
upsampleFactor = 10;
plotProperties = ["N", "P", "K", "OC", "Clay", "Sand", "Silt", "Moisture"];
ySlice = 0;
zSlice = []; % [] uses shallowest cell-center layer

krigingNeighbors = 32;
krigingRangeMeters = [12, 12, 1.5]; % [x y z] correlation ranges
krigingNuggetFraction = 1e-6;
predictionChunkSize = 750;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(scriptDir, 'scr')));

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
propertyFile = fullfile(outDir, 'staged_property_inversion_output.mat');
if ~isfile(propertyFile)
    error('Missing staged property output: %s', propertyFile);
end

S = load(propertyFile);
mesh = S.mesh;
config = S.config;
fields = localSelectFields(S, fieldSource);

xCenters = 0.5 * (mesh.x(1:end-1) + mesh.x(2:end));
yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
if isempty(zSlice)
    zSlice = zCenters(1);
end

xFine = linspace(xCenters(1), xCenters(end), numel(xCenters) * upsampleFactor);
yFine = linspace(yCenters(1), yCenters(end), numel(yCenters) * upsampleFactor);
zFine = linspace(zCenters(1), zCenters(end), numel(zCenters) * upsampleFactor);

[Xxz, Zxz] = meshgrid(xFine, zFine);
Yxz = ySlice * ones(size(Xxz));
targetXZ = [Xxz(:), Yxz(:), Zxz(:)];

[Xxy, Yxy] = meshgrid(xFine, yFine);
Zxy = zSlice * ones(size(Xxy));
targetXY = [Xxy(:), Yxy(:), Zxy(:)];

sourceXYZ = mesh.elementCenters;
kriged = struct();
kriged.xFine = xFine;
kriged.yFine = yFine;
kriged.zFine = zFine;
kriged.ySlice = ySlice;
kriged.zSlice = zSlice;
kriged.fieldSource = fieldSource;
kriged.upsampleFactor = upsampleFactor;
kriged.krigingNeighbors = krigingNeighbors;
kriged.krigingRangeMeters = krigingRangeMeters;

fprintf('\nKriging staged soil-property results\n');
fprintf('  Source: %s\n', propertyFile);
fprintf('  Field source: %s\n', fieldSource);
fprintf('  Source cells: %d\n', mesh.nElements);
fprintf('  Upsample factor: %d\n', upsampleFactor);
fprintf('  XZ target points: %d, XY target points: %d\n', size(targetXZ, 1), size(targetXY, 1));
fprintf('  Neighbors: %d, ranges [x y z] = [%.4g %.4g %.4g] m\n\n', ...
    krigingNeighbors, krigingRangeMeters);

for i = 1:numel(plotProperties)
    prop = char(plotProperties(i));
    if ~isfield(fields, prop)
        warning('Skipping %s: field is not present.', prop);
        continue;
    end

    sourceVals = fields.(prop)(:);
    fprintf('  Kriging %s ... ', prop);
    valsXZ = localOrdinaryKriging(sourceXYZ, sourceVals, targetXZ, ...
        krigingNeighbors, krigingRangeMeters, krigingNuggetFraction, predictionChunkSize);
    valsXY = localOrdinaryKriging(sourceXYZ, sourceVals, targetXY, ...
        krigingNeighbors, krigingRangeMeters, krigingNuggetFraction, predictionChunkSize);
    kriged.(prop).XZ = reshape(valsXZ, size(Xxz));
    kriged.(prop).XY = reshape(valsXY, size(Xxy));
    fprintf('done. XZ range [%.4g %.4g], XY range [%.4g %.4g]\n', ...
        min(valsXZ), max(valsXZ), min(valsXY), max(valsXY));
end

figXZ = localPlotKrigedGroup(kriged, plotProperties, "XZ", ...
    sprintf('Kriged soil properties, XZ slice (y ~= %.3g m)', ySlice), xFine, zFine);
figXY = localPlotKrigedGroup(kriged, plotProperties, "XY", ...
    sprintf('Kriged soil properties, XY slice (z ~= %.3g m)', zSlice), xFine, yFine);
figNCompare = localPlotOriginalVsKrigedN(mesh, fields.N, kriged.N, xFine, yFine, zFine, ySlice, zSlice);

if saveArg
    outMat = fullfile(outDir, 'staged_property_kriged_slices.mat');
    save(outMat, 'kriged', 'fieldSource', 'plotProperties', 'sourceXYZ', ...
        'krigingNeighbors', 'krigingRangeMeters', 'krigingNuggetFraction', '-v7.3');
    exportgraphics(figXZ, fullfile(outDir, 'ert_staged_property_kriged_xz_demo.png'), 'Resolution', 220);
    exportgraphics(figXY, fullfile(outDir, 'ert_staged_property_kriged_xy_demo.png'), 'Resolution', 220);
    exportgraphics(figNCompare, fullfile(outDir, 'ert_staged_property_kriged_N_comparison_demo.png'), 'Resolution', 220);
    fprintf('\nKriged outputs written to %s\n', outDir);
end

function fields = localSelectFields(S, fieldSource)
    switch lower(string(fieldSource))
        case "current"
            fields = S.fieldsProp;
        case "actual"
            if isfield(S, 'targetFieldsTrue') && ~isempty(S.targetFieldsTrue)
                fields = S.targetFieldsTrue;
            elseif isfield(S, 'fieldsStart')
                fields = S.fieldsStart;
            else
                error('No actual/target property fields found in output file.');
            end
        case "background"
            fields = S.fieldsPropertyBg;
        case "start"
            fields = S.fieldsStart;
        otherwise
            error('Unsupported fieldSource: %s', fieldSource);
    end
end

function predicted = localOrdinaryKriging(sourceXYZ, sourceVals, targetXYZ, nNeighbors, rangeMeters, nuggetFraction, chunkSize)
    sourceXYZ = double(sourceXYZ);
    sourceVals = double(sourceVals(:));
    targetXYZ = double(targetXYZ);
    rangeMeters = max(double(rangeMeters(:))', eps);
    nSource = size(sourceXYZ, 1);
    nNeighbors = min(nNeighbors, nSource);
    predicted = nan(size(targetXYZ, 1), 1);

    finiteMask = all(isfinite(sourceXYZ), 2) & isfinite(sourceVals);
    sourceXYZ = sourceXYZ(finiteMask, :);
    sourceVals = sourceVals(finiteMask);
    if isempty(sourceVals)
        return;
    end
    if max(sourceVals) == min(sourceVals)
        predicted(:) = sourceVals(1);
        return;
    end

    sill = var(sourceVals, 1);
    if ~(isfinite(sill) && sill > 0)
        sill = 1;
    end
    nugget = max(nuggetFraction * sill, eps);
    nTarget = size(targetXYZ, 1);

    scaledSource = sourceXYZ ./ rangeMeters;
    scaledTarget = targetXYZ ./ rangeMeters;
    for first = 1:chunkSize:nTarget
        last = min(first + chunkSize - 1, nTarget);
        for it = first:last
            d2 = sum((scaledSource - scaledTarget(it, :)).^2, 2);
            [~, order] = mink(d2, nNeighbors);
            xyzLocal = scaledSource(order, :);
            valsLocal = sourceVals(order);

            D = localPairwiseDistance(xyzLocal, xyzLocal);
            C = sill * exp(-D);
            C(1:size(C, 1)+1:end) = C(1:size(C, 1)+1:end) + nugget;

            dTarget = sqrt(sum((xyzLocal - scaledTarget(it, :)).^2, 2));
            c = sill * exp(-dTarget);

            A = [C, ones(nNeighbors, 1); ones(1, nNeighbors), 0];
            b = [c; 1];
            sol = A \ b;
            w = sol(1:nNeighbors);
            predicted(it) = sum(w .* valsLocal);
        end
    end
end

function D = localPairwiseDistance(A, B)
    D2 = max(sum(A.^2, 2) + sum(B.^2, 2)' - 2 * (A * B'), 0);
    D = sqrt(D2);
end

function fig = localPlotKrigedGroup(kriged, propNames, sliceName, figTitle, xVals, yVals)
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.04, 0.08, 0.92, 0.72]);
    tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, figTitle);
    for i = 1:numel(propNames)
        prop = char(propNames(i));
        nexttile(i);
        if ~isfield(kriged, prop)
            axis off;
            title(sprintf('%s missing', prop));
            continue;
        end
        values = kriged.(prop).(sliceName);
        imagesc(xVals, yVals, values);
        set(gca, 'YDir', 'normal');
        if sliceName == "XZ"
            set(gca, 'YDir', 'reverse');
            ylabel('z (m)');
        else
            ylabel('y (m)');
        end
        axis tight;
        colormap(gca, jet);
        colorbar;
        title(prop);
        xlabel('x (m)');
    end
end

function fig = localPlotOriginalVsKrigedN(mesh, nOriginal, nKriged, xFine, yFine, zFine, ySlice, zSlice)
    nXZ = localExtractXZSlice(mesh, nOriginal, ySlice);
    nXY = localExtractXYSlice(mesh, nOriginal, zSlice);
    climVals = localSharedClim([nXZ(:); nXY(:)], [nKriged.XZ(:); nKriged.XY(:)]);

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08, 0.12, 0.80, 0.68]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('N original vs kriged, field resolution x%d', ...
        round(numel(xFine) / max(1, numel(mesh.x) - 1))));

    nexttile(1);
    plot_rectilinear_cells(mesh.x, mesh.z, nXZ);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('Original N XZ (y ~= %.3g m)', ySlice));
    xlabel('x (m)'); ylabel('z (m)');

    nexttile(2);
    imagesc(xFine, zFine, nKriged.XZ);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('Kriged N XZ (y ~= %.3g m)', ySlice));
    xlabel('x (m)'); ylabel('z (m)');

    nexttile(3);
    plot_rectilinear_cells(mesh.x, mesh.y, nXY);
    set(gca, 'YDir', 'normal'); axis tight equal;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('Original N XY (z ~= %.3g m)', zSlice));
    xlabel('x (m)'); ylabel('y (m)');

    nexttile(4);
    imagesc(xFine, yFine, nKriged.XY);
    set(gca, 'YDir', 'normal'); axis tight equal;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('Kriged N XY (z ~= %.3g m)', zSlice));
    xlabel('x (m)'); ylabel('y (m)');
end

function vals2 = localExtractXZSlice(mesh, values, yTarget)
    vals3 = reshape(values, mesh.gridSize - 1);
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - yTarget));
    vals2 = squeeze(vals3(:, iy, :))';
end

function vals2 = localExtractXYSlice(mesh, values, zTarget)
    vals3 = reshape(values, mesh.gridSize - 1);
    zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
    [~, iz] = min(abs(zCenters - zTarget));
    vals2 = vals3(:, :, iz)';
end

function climVals = localSharedClim(varargin)
    vals = [];
    for i = 1:nargin
        vals = [vals; varargin{i}(:)]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
    if isempty(vals)
        climVals = [0, 1];
        return;
    end
    climVals = [min(vals), max(vals)];
    if climVals(1) == climVals(2)
        pad = max(abs(climVals(1)), 1) * 0.05;
        climVals = climVals + [-pad, pad];
    end
end
