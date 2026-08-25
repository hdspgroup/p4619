function fig = plot_ec_only_inversion_demo_3d(mesh, sigmaTrue, sigmaStart, sigmaCurrent, dObs, dPred, iterationLog, targetMisfit, summary, varargin)
%PLOT_EC_ONLY_INVERSION_DEMO_3D Slice-based summary figure for 3D EC inversion.

    p = inputParser;
    addParameter(p, 'ForwardMesh', []);
    addParameter(p, 'SigmaTrueForward', []);
    addParameter(p, 'SigmaTrueForwardNodes', []);
    addParameter(p, 'TrueDisplaySlices', []);
    addParameter(p, 'Survey', []);
    addParameter(p, 'TolObjective', []);
    parse(p, varargin{:});
    forwardMesh = p.Results.ForwardMesh;
    sigmaTrueForward = p.Results.SigmaTrueForward;
    sigmaTrueForwardNodes = p.Results.SigmaTrueForwardNodes;
    trueDisplaySlices = p.Results.TrueDisplaySlices;
    survey = p.Results.Survey;
    tolObjective = p.Results.TolObjective;
    cellSensitivity = [];
    if isfield(summary, 'cellSensitivity') && ~isempty(summary.cellSensitivity)
        cellSensitivity = summary.cellSensitivity(:);
    end

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.03, 0.08, 0.94, 0.78]);
    tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'EC-only inversion demo (3D)');

    ySlice = summary.trueModelMeta.sliceY;
    zSlice = summary.trueModelMeta.sliceZ;

    if ~isempty(trueDisplaySlices)
        xyActualVals = trueDisplaySlices.xy.v;
        xzActualVals = trueDisplaySlices.xz.v;
    elseif ~isempty(forwardMesh) && ~isempty(sigmaTrueForward)
        xyActualVals = localExtractXYSliceAny(forwardMesh, sigmaTrueForward, zSlice);
        xzActualVals = localExtractXZSliceAny(forwardMesh, sigmaTrueForward, ySlice);
    else
        xyActualVals = localExtractXYSlice(mesh, sigmaTrue, zSlice);
        xzActualVals = localExtractXZSlice(mesh, sigmaTrue, ySlice);
    end
    xyCurrentVals = localExtractXYSlice(mesh, sigmaCurrent, zSlice);
    xzCurrentVals = localExtractXZSlice(mesh, sigmaCurrent, ySlice);
    xyRange = localValueRange([xyActualVals(:); xyCurrentVals(:)]);
    xzRange = localValueRange([xzActualVals(:); xzCurrentVals(:)]);

    % Column 1: actual
    if ~isempty(trueDisplaySlices)
        nexttile(1); plotTrueDisplayXY(trueDisplaySlices.xy, 'EC actual XY (analytic slice)', xyRange, survey);
        nexttile(5); plotTrueDisplayXZ(trueDisplaySlices.xz, 'EC actual XZ (analytic slice)', xzRange);
    elseif ~isempty(forwardMesh) && ~isempty(sigmaTrueForward)
        nexttile(1); plotActualXYSlice(forwardMesh, sigmaTrueForward, sigmaTrueForwardNodes, zSlice, 'EC actual XY (fine, nodes)', xyRange, survey);
        nexttile(5); plotActualXZSlice(forwardMesh, sigmaTrueForward, ySlice, ...
            'EC actual XZ (fine, anomaly cross-section, nodes)', xzRange, sigmaTrueForwardNodes);
    else
        nexttile(1); plotXYSlice(mesh, sigmaTrue, zSlice, 'EC actual XY', xyRange, survey);
        nexttile(5); plotXZSlice(mesh, sigmaTrue, ySlice, ...
            'EC actual XZ (anomaly cross-section)', xzRange);
    end

    % Column 2: current inversion
    nexttile(2); plotXYSlice(mesh, sigmaCurrent, zSlice, 'EC current XY (coarse)', xyRange, [], cellSensitivity);
    nexttile(6); plotXZSlice(mesh, sigmaCurrent, ySlice, ...
        'EC current XZ (coarse, anomaly cross-section)', xzRange, cellSensitivity);

    % Column 3: residuals
    relResidual = (sigmaCurrent - sigmaTrue) ./ max(sigmaTrue, realmin);
    nexttile(3); plotCenteredResidualXYSlice(mesh, relResidual, zSlice, ...
        'EC residual XY ((current-actual)/actual)', cellSensitivity);
    nexttile(7); plotCenteredResidualXZSlice(mesh, relResidual, ySlice, ...
        'EC residual XZ ((current-actual)/actual)', cellSensitivity);

    % Column 4: diagnostics
    nexttile(4); plotVoltageScatter(dObs, dPred);
    nexttile(8); plotMisfitAndObjective(iterationLog, targetMisfit, tolObjective, summary);
end

function plotXYSlice(mesh, values, zSlice, ttl, climVals, survey, cellSensitivity)
    vals3 = reshape(values, mesh.gridSize - 1);
    iz = localNearestCellIndex(mesh.z, zSlice);
    cellVals = squeeze(vals3(:, :, iz))';
    alphaSlice = [];
    if nargin >= 7 && ~isempty(cellSensitivity)
        alphaSlice = localAlphaSliceXY(mesh, cellSensitivity, iz);
    end
    if isempty(alphaSlice)
        plot_rectilinear_cells(mesh.x, mesh.y, cellVals);
    else
        plot_rectilinear_cells(mesh.x, mesh.y, cellVals, 'AlphaData', alphaSlice);
    end
    axis equal tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 5 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (z ~= %.1f m)', ttl, localCellCenter(mesh.z, iz)));
    xlabel('x (m)');
    ylabel('y (m)');
    if nargin >= 6 && ~isempty(survey) && isfield(survey, 'electrodes')
        hold on;
        plot(survey.electrodes(:, 1), survey.electrodes(:, 2), 'ko', ...
            'MarkerSize', 4.5, 'LineWidth', 1.0, 'MarkerFaceColor', 'none');
        hold off;
    end
end

function plotActualXYSlice(mesh, values, nodalValues, zSlice, ttl, climVals, survey)
    if isfield(mesh, 'gridSize')
        plotXYSlice(mesh, values, zSlice, ttl, climVals, survey);
        return;
    end

    [xNodes, yNodes, vNodes, zUsed] = localExtractUnstructuredXYNodes(mesh, values, zSlice, nodalValues);
    tri = delaunay(xNodes, yNodes);
    trisurf(tri, xNodes, yNodes, zeros(size(xNodes)), vNodes, ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    axis equal tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 5 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (z ~= %.1f m)', ttl, zUsed));
    xlabel('x (m)');
    ylabel('y (m)');
    if nargin >= 6 && ~isempty(survey) && isfield(survey, 'electrodes')
        hold on;
        plot(survey.electrodes(:, 1), survey.electrodes(:, 2), 'ko', ...
            'MarkerSize', 4.5, 'LineWidth', 1.0, 'MarkerFaceColor', 'none');
        hold off;
    end
end

function plotTrueDisplayXY(slice, ttl, climVals, survey)
    surface(slice.x, slice.y, zeros(size(slice.v)), slice.v, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    axis equal tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 3 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (z ~= %.1f m)', ttl, slice.z));
    xlabel('x (m)');
    ylabel('y (m)');
    if nargin >= 4 && ~isempty(survey) && isfield(survey, 'electrodes')
        hold on;
        plot(survey.electrodes(:, 1), survey.electrodes(:, 2), 'ko', ...
            'MarkerSize', 4.5, 'LineWidth', 1.0, 'MarkerFaceColor', 'none');
        hold off;
    end
end

function plotXZSlice(mesh, values, ySlice, ttl, climVals, cellSensitivity)
    vals3 = reshape(values, mesh.gridSize - 1);
    iy = localNearestCellIndex(mesh.y, ySlice);
    cellVals = squeeze(vals3(:, iy, :))';
    alphaSlice = [];
    if nargin >= 6 && ~isempty(cellSensitivity)
        alphaSlice = localAlphaSliceXZ(mesh, cellSensitivity, iy);
    end
    if isempty(alphaSlice)
        plot_rectilinear_cells(mesh.x, mesh.z, cellVals);
    else
        plot_rectilinear_cells(mesh.x, mesh.z, cellVals, 'AlphaData', alphaSlice);
    end
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 5 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (y ~= %.1f m)', ttl, localCellCenter(mesh.y, iy)));
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotActualXZSlice(mesh, values, ySlice, ttl, climVals, nodalValues)
    if isfield(mesh, 'gridSize')
        plotXZSlice(mesh, values, ySlice, ttl, climVals);
        return;
    end

    [xNodes, zNodes, vNodes, yUsed] = localExtractUnstructuredXZNodes(mesh, values, ySlice, nodalValues);
    tri = delaunay(xNodes, zNodes);
    trisurf(tri, xNodes, zNodes, zeros(size(xNodes)), vNodes, ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 5 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (y ~= %.1f m)', ttl, yUsed));
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotTrueDisplayXZ(slice, ttl, climVals)
    surface(slice.x, slice.z, zeros(size(slice.v)), slice.v, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 3 && ~isempty(climVals)
        clim(climVals);
    end
    title(sprintf('%s (y ~= %.1f m)', ttl, slice.y));
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotCenteredResidualXYSlice(mesh, values, zSlice, ttl, cellSensitivity)
    vals3 = reshape(values, mesh.gridSize - 1);
    iz = localNearestCellIndex(mesh.z, zSlice);
    sliceVals = squeeze(vals3(:, :, iz))';
    alphaSlice = [];
    if nargin >= 5 && ~isempty(cellSensitivity)
        alphaSlice = localAlphaSliceXY(mesh, cellSensitivity, iz);
    end
    if isempty(alphaSlice)
        plot_rectilinear_cells(mesh.x, mesh.y, sliceVals);
    else
        plot_rectilinear_cells(mesh.x, mesh.y, sliceVals, 'AlphaData', alphaSlice);
    end
    axis equal tight;
    colormap(gca, jet);
    cmax = localResidualRange(sliceVals);
    clim([-cmax, cmax]);
    cb = colorbar;
    ylabel(cb, '(current - actual) / actual');
    title(sprintf('%s (z ~= %.1f m)', ttl, localCellCenter(mesh.z, iz)));
    xlabel('x (m)');
    ylabel('y (m)');
end

function plotCenteredResidualXZSlice(mesh, values, ySlice, ttl, cellSensitivity)
    vals3 = reshape(values, mesh.gridSize - 1);
    iy = localNearestCellIndex(mesh.y, ySlice);
    sliceVals = squeeze(vals3(:, iy, :))';
    alphaSlice = [];
    if nargin >= 5 && ~isempty(cellSensitivity)
        alphaSlice = localAlphaSliceXZ(mesh, cellSensitivity, iy);
    end
    if isempty(alphaSlice)
        plot_rectilinear_cells(mesh.x, mesh.z, sliceVals);
    else
        plot_rectilinear_cells(mesh.x, mesh.z, sliceVals, 'AlphaData', alphaSlice);
    end
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap(gca, jet);
    cmax = localResidualRange(sliceVals);
    clim([-cmax, cmax]);
    cb = colorbar;
    ylabel(cb, '(current - actual) / actual');
    title(sprintf('%s (y ~= %.1f m)', ttl, localCellCenter(mesh.y, iy)));
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotVoltageScatter(dObs, dPred)
    hold on;
    scatter(dObs, dPred, 30, 'filled', 'MarkerFaceColor', [0.15 0.35 0.85], ...
        'MarkerFaceAlpha', 0.7, 'DisplayName', 'measurements');
    lims = [min([dObs(:); dPred(:)]), max([dObs(:); dPred(:)])];
    if ~all(isfinite(lims)) || diff(lims) <= 0
        lims = [-1, 1];
    end
    plot(lims, lims, 'k--', 'LineWidth', 1.2, 'DisplayName', '1:1');
    xlim(lims);
    ylim(lims);
    axis square;
    grid on;
    xlabel('observed voltage (V)');
    ylabel('current voltage (V)');
    title('Voltage fit');
    legend('Location', 'best');
end

function plotMisfitAndObjective(iterationLog, targetMisfit, tolObjective, summary)
    it = [iterationLog.iter];
    misfitVals = max([iterationLog.normalizedDataMisfit], realmin);
    objVals = max([iterationLog.objective], realmin);
    objValsNorm = objVals / max(objVals(1), realmin);

    hold on;
    h1 = semilogy(it, misfitVals, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', ...
        'DisplayName', 'data misfit');
    h2 = semilogy(it, objValsNorm, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'total objective / initial');
    h3 = yline(max(targetMisfit, realmin), 'b--', 'LineWidth', 1.3, ...
        'DisplayName', 'target misfit');
    if ~isempty(tolObjective)
        h4 = yline(max(tolObjective, realmin), 'r--', 'LineWidth', 1.3, ...
            'DisplayName', 'objective tol');
        legend([h1, h2, h3, h4], 'Location', 'best');
    else
        legend([h1, h2, h3], 'Location', 'best');
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('iteration');
    ylabel('log scale');
    title(sprintf('Convergence (%s)', summary.exitReason));
end

function idx = localNearestCellIndex(axisEdges, targetVal)
    centers = 0.5 * (axisEdges(1:end-1) + axisEdges(2:end));
    [~, idx] = min(abs(centers - targetVal));
end

function sliceVals = localExtractXYSliceAny(mesh, values, zSlice)
    if isfield(mesh, 'gridSize')
        sliceVals = localExtractXYSlice(mesh, values, zSlice);
    else
        [~, ~, sliceVals] = localSampleUnstructuredXY(mesh, values, zSlice);
    end
end

function sliceVals = localExtractXZSliceAny(mesh, values, ySlice)
    if isfield(mesh, 'gridSize')
        sliceVals = localExtractXZSlice(mesh, values, ySlice);
    else
        [~, ~, sliceVals] = localSampleUnstructuredXZ(mesh, values, ySlice);
    end
end

function sliceVals = localExtractXYSlice(mesh, values, zSlice)
    vals3 = reshape(values, mesh.gridSize - 1);
    iz = localNearestCellIndex(mesh.z, zSlice);
    sliceVals = squeeze(vals3(:, :, iz))';
end

function sliceVals = localExtractXZSlice(mesh, values, ySlice)
    vals3 = reshape(values, mesh.gridSize - 1);
    iy = localNearestCellIndex(mesh.y, ySlice);
    sliceVals = squeeze(vals3(:, iy, :))';
end

function vr = localValueRange(vals)
    vals = vals(isfinite(vals));
    if isempty(vals)
        vr = [0, 1];
        return;
    end

    vmin = min(vals);
    vmax = max(vals);
    if ~(isfinite(vmin) && isfinite(vmax)) || vmax <= vmin
        pad = max(abs(vmin), 1) * 0.05;
        vr = [vmin - pad, vmin + pad];
        return;
    end

    vr = [vmin, vmax];
end

function [Xq, Yq, Vq] = localSampleUnstructuredXY(mesh, values, zSlice)
    c = mesh.elementCenters;
    x = c(:, 1); y = c(:, 2); z = c(:, 3);
    nx = 140; ny = 140;
    xq = linspace(min(x), max(x), nx);
    yq = linspace(min(y), max(y), ny);
    [Xq, Yq] = meshgrid(xq, yq);
    F = scatteredInterpolant(x, y, z, values(:), 'natural', 'nearest');
    Vq = F(Xq, Yq, zSlice * ones(size(Xq)));
end

function [Xq, Zq, Vq] = localSampleUnstructuredXZ(mesh, values, ySlice)
    c = mesh.elementCenters;
    x = c(:, 1); y = c(:, 2); z = c(:, 3);
    nx = 140; nz = 120;
    xq = linspace(min(x), max(x), nx);
    zq = linspace(min(z), max(z), nz);
    [Xq, Zq] = meshgrid(xq, zq);
    F = scatteredInterpolant(x, y, z, values(:), 'natural', 'nearest');
    Vq = F(Xq, ySlice * ones(size(Xq)), Zq);
end

function [xNodes, yNodes, vNodes, zUsed] = localExtractUnstructuredXYNodes(mesh, values, zSlice, nodalValues)
    if nargin < 4 || isempty(nodalValues)
        nodalVals = localElementToNodeValues(mesh, values);
    else
        nodalVals = nodalValues(:);
    end
    zAll = mesh.nodes(:, 3);
    zUnique = unique(round(zAll, 8));
    [~, iz] = min(abs(zUnique - zSlice));
    zUsed = zUnique(iz);
    mask = abs(zAll - zUsed) < 1e-8;
    xNodes = mesh.nodes(mask, 1);
    yNodes = mesh.nodes(mask, 2);
    vNodes = nodalVals(mask);
end

function [xNodes, zNodes, vNodes, yUsed] = localExtractUnstructuredXZNodes(mesh, values, ySlice, nodalValues)
    if nargin < 4 || isempty(nodalValues)
        nodalVals = localElementToNodeValues(mesh, values);
    else
        nodalVals = nodalValues(:);
    end
    yAll = mesh.nodes(:, 2);
    yUnique = unique(sort(round(yAll, 8)));
    if numel(yUnique) > 1
        dy = diff(yUnique);
        dy = dy(dy > 1e-8);
        if isempty(dy)
            bandTol = 1e-8;
        else
            bandTol = max(0.5 * min(dy), 0.75);
        end
    else
        bandTol = 1e-8;
    end
    [~, iy] = min(abs(yUnique - ySlice));
    yUsed = yUnique(iy);
    mask = abs(yAll - yUsed) <= bandTol;
    xNodes = mesh.nodes(mask, 1);
    zNodes = mesh.nodes(mask, 3);
    vNodes = nodalVals(mask);
end

function nodalVals = localElementToNodeValues(mesh, elemVals)
    nNodes = size(mesh.nodes, 1);
    nElem = size(mesh.elements, 1);
    rows = mesh.elements(:);
    cols = repelem((1:nElem)', size(mesh.elements, 2), 1);
    A = sparse(rows, cols, 1, nNodes, nElem);
    counts = full(sum(A, 2));
    counts(counts == 0) = 1;
    nodalVals = (A * elemVals(:)) ./ counts;
end

function centerVal = localCellCenter(axisEdges, idx)
    centers = 0.5 * (axisEdges(1:end-1) + axisEdges(2:end));
    centerVal = centers(idx);
end

function cmax = localResidualRange(sliceVals)
    cmax = max(abs(sliceVals(:)));
    if ~isfinite(cmax) || cmax == 0
        cmax = 1;
    end
end

function alphaSlice = localAlphaSliceXY(mesh, cellSensitivity, iz)
    sens3 = reshape(cellSensitivity, mesh.gridSize - 1);
    alphaSlice = squeeze(sens3(:, :, iz))';
    alphaSlice = localNormalizeAlpha(alphaSlice);
end

function alphaSlice = localAlphaSliceXZ(mesh, cellSensitivity, iy)
    sens3 = reshape(cellSensitivity, mesh.gridSize - 1);
    alphaSlice = squeeze(sens3(:, iy, :))';
    alphaSlice = localNormalizeAlpha(alphaSlice);
end

function alphaVals = localNormalizeAlpha(alphaVals)
    alphaVals(~isfinite(alphaVals)) = 0;
    if all(alphaVals(:) == 0)
        alphaVals = ones(size(alphaVals));
        return;
    end

    alphaVals = alphaVals / max(alphaVals(:));
    alphaVals = sqrt(alphaVals);
    alphaVals = 0.15 + 0.85 * alphaVals;
end
