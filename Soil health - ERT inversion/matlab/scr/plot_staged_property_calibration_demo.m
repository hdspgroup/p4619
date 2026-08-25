function figs = plot_staged_property_calibration_demo(mesh, config, fieldsTrue, fieldsCurrent, sigmaTarget, sigmaCurrent, ...
    iterationLog, summary, refFields, aux, ecModel, truthSettings, activeMask, maxScaleFactor, residualOverviewScaleMode, residualColorLimit)
%PLOT_STAGED_PROPERTY_CALIBRATION_DEMO Plot staged soil-property calibration.

    if nargin < 13 || isempty(activeMask)
        activeMask = true(mesh.nElements, 1);
    end
    if nargin < 14 || isempty(maxScaleFactor)
        maxScaleFactor = 25;
    end
    if nargin < 15 || strlength(string(residualOverviewScaleMode)) == 0
        residualOverviewScaleMode = "shared";
    end
    if nargin < 16
        residualColorLimit = [];
    end

    ySlice = 0;
    xDisplay = linspace(mesh.x(1), mesh.x(end), 241);
    zDisplay = linspace(mesh.z(1), mesh.z(end), 161);

    trueCellFields = evaluate_soil_true_model_at_points(mesh.elementCenters, refFields, truthSettings);
    sigmaBestXZ = localExtractXZSlice(mesh, sigmaTarget, ySlice);
    sigmaPropXZ = localExtractXZSlice(mesh, sigmaCurrent, ySlice);

    figs = struct();
    figs.npkOc = localPlotPropertyGroup(mesh, trueCellFields, fieldsCurrent, ...
        xDisplay, zDisplay, ["N", "P", "K", "OC"], "N/P/K/OC property calibration (XZ, y ~= 0 m)", activeMask, maxScaleFactor, residualColorLimit);
    figs.textureMoisture = localPlotPropertyGroup(mesh, trueCellFields, fieldsCurrent, ...
        xDisplay, zDisplay, ["Clay", "Sand", "Silt", "Moisture"], "Texture/moisture property calibration (XZ, y ~= 0 m)", activeMask, maxScaleFactor, residualColorLimit);
    figs.residualOverview = localPlotResidualOverview(mesh, trueCellFields, fieldsCurrent, ...
        xDisplay, zDisplay, ["N", "P", "K", "OC", "Clay", "Sand", "Silt", "Moisture"], activeMask, residualOverviewScaleMode, residualColorLimit);
    figs.diagnostics = localPlotDiagnostics(mesh, sigmaBestXZ, sigmaPropXZ, ...
        xDisplay, zDisplay, sigmaTarget, sigmaCurrent, iterationLog, summary, activeMask, maxScaleFactor);
end

function fig = localPlotPropertyGroup(mesh, trueCellFields, fieldsCurrent, xDisplay, zDisplay, propNames, figTitle, activeMask, maxScaleFactor, residualColorLimit)
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.04, 0.08, 0.92, 0.78]);
    tl = tiledlayout(3, numel(propNames), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, figTitle);
    activeXZ = localExtractXZSlice(mesh, double(activeMask(:)), 0) > 0.5;
    alphaXZ = localActiveAlpha(activeXZ);

    for i = 1:numel(propNames)
        prop = char(propNames(i));
        actualGrid = localExtractXZSlice(mesh, trueCellFields.(prop), 0);
        currentGrid = localExtractXZSlice(mesh, fieldsCurrent.(prop), 0);
        residualGrid = localRelativeResidual(currentGrid, actualGrid);

        climVals = localAlignedLinearClim(actualGrid, currentGrid, maxScaleFactor);
        fprintf('  Plot color limits %-10s: [%.4g %.4g]\n', prop, climVals(1), climVals(2));
        [resDisp, resClim, resTicks, resLabels] = localSignedLogResidualData(residualGrid, residualColorLimit);

        nexttile(i);
        localPlotCurrentSlice(mesh.x, mesh.z, actualGrid, climVals, sprintf('%s actual', prop), [], [], alphaXZ, activeXZ);

        nexttile(i + numel(propNames));
        localPlotCurrentSlice(mesh.x, mesh.z, currentGrid, climVals, sprintf('%s current', prop), [], [], alphaXZ, activeXZ);

        nexttile(i + 2 * numel(propNames));
        localPlotResidualSlice(mesh.x, mesh.z, resDisp, resClim, ...
            sprintf('%s residual ((current-actual)/actual)', prop), resTicks, resLabels, alphaXZ, activeXZ);
    end
end

function fig = localPlotResidualOverview(mesh, trueCellFields, fieldsCurrent, xDisplay, zDisplay, propNames, activeMask, scaleMode, residualColorLimit)
    residuals = cell(numel(propNames), 1);
    allResiduals = [];
    axList = gobjects(numel(propNames), 1);
    activeXZ = localExtractXZSlice(mesh, double(activeMask(:)), 0) > 0.5;
    alphaXZ = localActiveAlpha(activeXZ);
    [fixedResidualLimit, ~] = localParseResidualColorLimit(residualColorLimit);
    for i = 1:numel(propNames)
        prop = char(propNames(i));
        actualGrid = localExtractXZSlice(mesh, trueCellFields.(prop), 0);
        currentGrid = localExtractXZSlice(mesh, fieldsCurrent.(prop), 0);
        residuals{i} = localRelativeResidual(currentGrid, actualGrid);
        allResiduals = [allResiduals; residuals{i}(:)]; %#ok<AGROW>
    end

    scaleMode = lower(string(scaleMode));
    if scaleMode == "shared"
        [~, sharedClimVals, sharedTicks, sharedLabels, sharedScale] = localSignedLogResidualData(allResiduals, residualColorLimit);
        if ~fixedResidualLimit
            [sharedTicks, sharedLabels] = localSignedLogOverviewTicks(sharedScale, sharedClimVals);
        end
    elseif scaleMode ~= "perproperty"
        error('Unsupported residual overview scale mode: %s', char(scaleMode));
    end

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.04, 0.08, 0.92, 0.70]);
    tl = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Soil-property residual overview (XZ, y ~= 0 m, %s scale)', char(scaleMode)));

    for i = 1:numel(propNames)
        prop = char(propNames(i));
        if scaleMode == "shared"
            if fixedResidualLimit
                resDisp = residuals{i};
            else
                resDisp = localSignedLogTransform(residuals{i}, sharedScale);
            end
            climVals = sharedClimVals;
            ticks = sharedTicks;
            labels = sharedLabels;
        else
            [resDisp, climVals, ticks, labels, localScale] = localSignedLogResidualData(residuals{i}, residualColorLimit);
            if ~fixedResidualLimit
                [ticks, labels] = localSignedLogOverviewTicks(localScale, climVals);
            end
        end
        fprintf('  Residual overview limits %-10s: [%.4g %.4g] (%s)\n', prop, climVals(1), climVals(2), char(scaleMode));
        axList(i) = nexttile(i);
        localPlotResidualSlice(mesh.x, mesh.z, resDisp, climVals, ...
            sprintf('%s residual', prop), ticks, labels, alphaXZ, activeXZ);
    end

    if scaleMode == "shared"
        for i = 1:numel(axList)
            if isgraphics(axList(i))
                colormap(axList(i), jet);
                clim(axList(i), sharedClimVals);
            end
        end
    end
end

function fig = localPlotDiagnostics(mesh, sigmaBestXZ, sigmaPropXZ, xDisplay, zDisplay, ...
        sigmaTarget, sigmaCurrent, iterationLog, summary, activeMask, maxScaleFactor)

    sigmaClim = localAlignedLinearClim(sigmaBestXZ, sigmaPropXZ, maxScaleFactor);
    fprintf('  Plot color limits %-10s: [%.4g %.4g]\n', 'EC', sigmaClim(1), sigmaClim(2));
    sigmaResidual = localRelativeResidual(sigmaPropXZ, sigmaBestXZ);
    [resSigmaDisp, resSigmaClim, resSigmaTicks, resSigmaLabels] = localSignedLogResidualData(sigmaResidual);
    activeXZ = localExtractXZSlice(mesh, double(activeMask(:)), 0) > 0.5;
    alphaXZ = localActiveAlpha(activeXZ);

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08, 0.12, 0.84, 0.68]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Staged soil-property calibration diagnostics');

    nexttile(1);
    localPlotCurrentSlice(mesh.x, mesh.z, sigmaBestXZ, sigmaClim, 'EC target (cell XZ)', [], [], alphaXZ, activeXZ);

    nexttile(2);
    localPlotCurrentSlice(mesh.x, mesh.z, sigmaPropXZ, sigmaClim, 'Property-derived EC current (XZ)', [], [], alphaXZ, activeXZ);

    nexttile(3);
    localPlotResidualSlice(mesh.x, mesh.z, resSigmaDisp, resSigmaClim, ...
        'EC residual ((propertyEC-targetEC)/targetEC)', resSigmaTicks, resSigmaLabels, alphaXZ, activeXZ);

    nexttile(4);
    localPlotEcScatter(sigmaTarget(activeMask), sigmaCurrent(activeMask));

    nexttile(5);
    localPlotMisfitHistory(iterationLog);

    nexttile(6);
    localPlotRunSummary(summary);
end

function localPlotActualSlice(xDisplay, zDisplay, values, climVals, ttl, cbTicks, cbLabels)
    imagesc(xDisplay, zDisplay, values);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks), set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels); end
    title(ttl); xlabel('x (m)'); ylabel('z (m)');
end

function localPlotCurrentSlice(xEdges, zEdges, cellValues, climVals, ttl, cbTicks, cbLabels, alphaData, activeMask2)
    if nargin < 8, alphaData = []; end
    if nargin < 9, activeMask2 = []; end
    plot_rectilinear_cells(xEdges, zEdges, cellValues, 'AlphaData', alphaData);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks), set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels); end
    localDrawActiveBoundary(xEdges, zEdges, activeMask2);
    title(ttl); xlabel('x (m)'); ylabel('z (m)');
end

function localPlotResidualSlice(xEdges, zEdges, cellValues, climVals, ttl, cbTicks, cbLabels, alphaData, activeMask2)
    if nargin < 8, alphaData = []; end
    if nargin < 9, activeMask2 = []; end
    plot_rectilinear_cells(xEdges, zEdges, cellValues, 'AlphaData', alphaData);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks), set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels); end
    localDrawActiveBoundary(xEdges, zEdges, activeMask2);
    title(ttl); xlabel('x (m)'); ylabel('z (m)');
end

function localPlotEcScatter(sigmaTarget, sigmaCurrent)
    targetLog = log10(max(sigmaTarget(:), realmin));
    currentLog = log10(max(sigmaCurrent(:), realmin));
    scatter(targetLog, currentLog, 18, 'filled', ...
        'MarkerFaceColor', [0.15 0.35 0.85], 'MarkerFaceAlpha', 0.55, 'DisplayName', 'cells');
    hold on;
    lims = [min([targetLog; currentLog]), max([targetLog; currentLog])];
    if ~all(isfinite(lims)) || diff(lims) <= 0, lims = [-1 1]; end
    pad = 0.05 * max(diff(lims), 1);
    lims = lims + [-pad pad];
    plot(lims, lims, 'k--', 'LineWidth', 1.2, 'DisplayName', '1:1');
    xlim(lims); ylim(lims); axis square; grid on;
    xlabel('log_{10} target EC'); ylabel('log_{10} property EC');
    title('Cell EC match');
    legend('Location', 'best');
end

function alphaData = localActiveAlpha(activeMask2)
    alphaData = 0.22 * ones(size(activeMask2));
    alphaData(activeMask2) = 1.0;
end

function localDrawActiveBoundary(xEdges, zEdges, activeMask2)
    if isempty(activeMask2)
        return;
    end
    hold on;
    boundaryVals = double(activeMask2);
    xCenters = 0.5 * (xEdges(1:end-1) + xEdges(2:end));
    zCenters = 0.5 * (zEdges(1:end-1) + zEdges(2:end));
    if numel(unique(boundaryVals(:))) > 1
        contour(xCenters, zCenters, boundaryVals, [0.5 0.5], ...
            'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

function localPlotMisfitHistory(iterationLog)
    if isempty(iterationLog)
        axis off; title('Misfit history'); text(0.1, 0.5, 'No property iterations logged');
        return;
    end
    it = [iterationLog.iter];
    semilogy(it, max([iterationLog.normalizedDataMisfit], realmin), 'b-o', ...
        'LineWidth', 1.4, 'MarkerFaceColor', 'b');
    grid on; xlabel('iteration'); ylabel('normalized EC misfit');
    title('Property EC convergence');
end

function localPlotRunSummary(summary)
    axis off;
    txt = sprintf(['Exit: %s\nIterations: %d\nObj before: %.3g\nObj after: %.3g\n' ...
        'Misfit after: %.3g\nPrior rows: %d\nPrior penalty: %.3g\nPrior EC penalty: %.3g\nStep length: %.3g'], ...
        summary.exitReason, summary.iterationsCompleted, summary.objectiveBefore, summary.objectiveAfter, ...
        summary.normalizedDataMisfitAfter, summary.priorRows, summary.priorPenaltyAfter, ...
        summary.priorEcPenaltyAfter, summary.lastStepLength);
    text(0.05, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 11, 'FontWeight', 'bold');
    title('Run summary');
end

function vals2 = localExtractXZSlice(mesh, values, yTarget)
    vals3 = reshape(values, mesh.gridSize - 1);
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - yTarget));
    vals2 = squeeze(vals3(:, iy, :))';
end

function sampled = localSampleActualToCoarse(mesh, actualGrid, xDisplay, zDisplay)
    xCenters = 0.5 * (mesh.x(1:end-1) + mesh.x(2:end));
    zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
    [Xc, Zc] = meshgrid(xCenters, zCenters);
    sampled = interp2(xDisplay, zDisplay, actualGrid, Xc, Zc, 'linear');
end

function sigma = localPredictSigmaFromFields(fields, config, aux, ecModel)
    nCells = numel(fields.N);
    cfg = localResizeConfig(config, nCells);
    auxEval = struct();
    auxNames = fieldnames(aux);
    for i = 1:numel(auxNames)
        val = aux.(auxNames{i});
        if isscalar(val)
            auxEval.(auxNames{i}) = repmat(val, nCells, 1);
        else
            auxEval.(auxNames{i}) = val(1) * ones(nCells, 1);
        end
    end
    m = pack_soil_fields(fields, cfg);
    sigma = predict_soil_ec_physics(m, cfg, auxEval, ecModel);
end

function cfg = localResizeConfig(config, nCells)
    cfg = config;
    cfg.nCells = nCells;
    cfg.totalSize = config.nProps * nCells;
    index = struct();
    for i = 1:config.nProps
        first = (i - 1) * nCells + 1;
        last = i * nCells;
        index.(char(config.names(i))) = first:last;
    end
    cfg.index = index;
end

function ratio = localRelativeResidual(currentVals, actualVals)
    denom = max(abs(actualVals), 1e-6);
    ratio = (currentVals - actualVals) ./ denom;
end

function climVals = localAlignedLinearClim(actualVals, currentVals, maxScaleFactor)
    actualFinite = actualVals(isfinite(actualVals));
    currentFinite = currentVals(isfinite(currentVals));
    if isempty(actualFinite) && isempty(currentFinite), climVals = [0, 1]; return; end
    if isempty(actualFinite), actualFinite = currentFinite; end
    if isempty(currentFinite), currentFinite = actualFinite; end
    combinedMin = min([actualFinite; currentFinite]); combinedMax = max([actualFinite; currentFinite]);
    if ~isfinite(maxScaleFactor) || maxScaleFactor <= 0
        lower = combinedMin;
        upper = combinedMax;
    else
        actualMin = min(actualFinite); actualMax = max(actualFinite); actualSpan = actualMax - actualMin;
        if ~(isfinite(actualSpan) && actualSpan > 0), actualSpan = max(abs(actualFinite)); end
        if ~(isfinite(actualSpan) && actualSpan > 0), actualSpan = 1; end
        allowedSpan = maxScaleFactor * actualSpan;
        lower = combinedMin; upper = combinedMax;
        if (upper - lower) > allowedSpan
            lower = actualMin; upper = lower + allowedSpan;
            if combinedMax > upper, upper = combinedMax; lower = upper - allowedSpan; end
        end
    end
    if lower == upper
        pad = max(abs(lower), 1) * 0.5;
        climVals = [lower - pad, upper + pad];
    else
        climVals = [lower, upper];
    end
end

function [dispVals, climVals, ticks, labels, scale] = localSignedLogResidualData(residualVals, residualColorLimit)
    if nargin < 2
        residualColorLimit = [];
    end
    [fixedLimit, fixedClim] = localParseResidualColorLimit(residualColorLimit);
    if fixedLimit
        dispVals = residualVals;
        climVals = fixedClim;
        ticks = localDirectResidualTicks(climVals);
        labels = arrayfun(@localFormatResidualTick, ticks, 'UniformOutput', false);
        scale = 1;
        return;
    end
    finiteVals = residualVals(isfinite(residualVals));
    if isempty(finiteVals)
        dispVals = residualVals; climVals = [-1, 1]; ticks = []; labels = {}; scale = 1; return;
    end
    absVals = abs(finiteVals);
    absVals = absVals(absVals > 0);
    if isempty(absVals)
        dispVals = zeros(size(residualVals)); climVals = [-1, 1];
        scale = 1;
        [ticks, labels] = localSignedLogOverviewTicks(scale, climVals);
        return;
    end
    scale = max(prctile(absVals, 25), 1e-6);
    dispVals = localSignedLogTransform(residualVals, scale);
    cmax = max(abs(dispVals(isfinite(dispVals))));
    if ~(isfinite(cmax) && cmax > 0), cmax = 1; end
    climVals = [-cmax, cmax];
    [ticks, labels] = localSignedLogOverviewTicks(scale, climVals);
end

function [isFixed, climVals] = localParseResidualColorLimit(residualColorLimit)
    isFixed = false;
    climVals = [];
    if isempty(residualColorLimit) || ~isnumeric(residualColorLimit)
        return;
    end
    vals = residualColorLimit(:)';
    vals = vals(isfinite(vals));
    if isscalar(vals) && vals > 0
        climVals = [-abs(vals), abs(vals)];
        isFixed = true;
    elseif numel(vals) == 2 && vals(1) < vals(2)
        climVals = vals;
        isFixed = true;
    end
end

function ticks = localDirectResidualTicks(climVals)
    lo = climVals(1);
    hi = climVals(2);
    if lo < 0 && hi > 0
        mag = max(abs(climVals));
        baseTicks = [-mag, -0.3 * mag, -0.1 * mag, 0, 0.1 * mag, 0.3 * mag, mag];
    else
        baseTicks = linspace(lo, hi, 7);
    end
    ticks = unique(baseTicks(baseTicks >= lo - eps & baseTicks <= hi + eps), 'stable');
    if numel(ticks) > 7
        ticks = ticks(round(linspace(1, numel(ticks), 7)));
    end
end

function [ticks, labels] = localSignedLogOverviewTicks(scale, climVals)
    cmax = max(abs(climVals(:)));
    if ~(isfinite(cmax) && cmax > 0)
        ticks = [-1 0 1];
        labels = {'-1','0','1'};
        return;
    end

    maxResidual = max(scale, 1e-12) * (10 ^ cmax - 1);
    if ~(isfinite(maxResidual) && maxResidual > 0)
        ticks = [0];
        labels = {'0'};
        return;
    end

    expMin = floor(log10(maxResidual)) - 4;
    expMax = ceil(log10(maxResidual));
    magCandidates = [];
    for e = expMin:expMax
        magCandidates = [magCandidates, 10 ^ e, 3 * 10 ^ e]; %#ok<AGROW>
    end
    magCandidates = unique(magCandidates(isfinite(magCandidates) & magCandidates > 0 & magCandidates <= maxResidual * (1 + 1e-10)));
    if isempty(magCandidates)
        magCandidates = maxResidual;
    end
    magCandidates = sort(magCandidates(:));

    % Seven total ticks means at most three positive magnitudes plus their
    % negative counterparts and zero.
    if numel(magCandidates) > 3
        magCandidates = magCandidates(end-2:end);
    end

    residualTicks = [-flipud(magCandidates); 0; magCandidates];
    ticks = localSignedLogTransform(residualTicks, scale);
    keep = ticks >= climVals(1) - 1e-10 & ticks <= climVals(2) + 1e-10;
    residualTicks = residualTicks(keep);
    ticks = ticks(keep);
    labels = arrayfun(@localFormatResidualTick, residualTicks, 'UniformOutput', false);
end

function label = localFormatResidualTick(x)
    if x == 0
        label = '0';
    elseif abs(x) >= 1
        label = sprintf('%.0f', x);
    elseif abs(x) >= 0.1
        label = sprintf('%.1f', x);
    elseif abs(x) >= 0.01
        label = sprintf('%.2f', x);
    else
        label = sprintf('%.3g', x);
    end
end

function dispVals = localSignedLogTransform(values, scale)
    dispVals = sign(values) .* log10(1 + abs(values) ./ max(scale, 1e-12));
end
