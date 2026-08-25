function fig = plot_property_inversion_demo(mesh, config, ...
    fieldsTrue, fieldsCurrent, sigmaTrue, sigmaCurrent, ...
    dObs, dPred, iterationLog, targetMisfit, summary, ...
    refFields, aux, ecModel, truthSettings)
%PLOT_PROPERTY_INVERSION_DEMO Property inversion summary in EC-only style.

    if nargin < 14
        truthSettings = [];
    end

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'normalized', ...
        'Position', [0.03, 0.08, 0.94, 0.82]);
    tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Property-decomposed ERT inversion demo');

    ySlice = 0;
    xDisplay = linspace(mesh.x(1), mesh.x(end), 241);
    zDisplay = linspace(mesh.z(1), mesh.z(end), 161);
    [Xg, Zg] = meshgrid(xDisplay, zDisplay);
    Yg = ySlice * ones(size(Xg));
    slicePoints = [Xg(:), Yg(:), Zg(:)];

    trueSliceFields = evaluate_soil_true_model_at_points(slicePoints, refFields, truthSettings);
    sigmaSlice = localPredictSigmaFromFields(trueSliceFields, config, aux, ecModel);

    actualEC = reshape(sigmaSlice, size(Xg));
    actualN = reshape(trueSliceFields.N, size(Xg));
    actualClay = reshape(trueSliceFields.Clay, size(Xg));

    currentEC = localExtractXZSlice(mesh, sigmaCurrent, ySlice);
    currentN = localExtractXZSlice(mesh, fieldsCurrent.N, ySlice);
    currentClay = localExtractXZSlice(mesh, fieldsCurrent.Clay, ySlice);

    residualEC = localRelativeResidual(currentEC, localSampleActualToCoarse(mesh, actualEC, xDisplay, zDisplay, ySlice));
    residualN = localRelativeResidual(currentN, localSampleActualToCoarse(mesh, actualN, xDisplay, zDisplay, ySlice));
    residualClay = localRelativeResidual(currentClay, localSampleActualToCoarse(mesh, actualClay, xDisplay, zDisplay, ySlice));
    [residualECDisp, residualECClim, residualECTicks, residualECLabels] = localSignedLogResidualData(residualEC);
    [residualNDisp, residualNClim, residualNTicks, residualNLabels] = localSignedLogResidualData(residualN);
    [residualClayDisp, residualClayClim, residualClayTicks, residualClayLabels] = localSignedLogResidualData(residualClay);

    maxScaleFactor = 5;
    [actualECDisp, currentECDisp, ecClim, ecTicks, ecLabels] = localLogDisplayData(actualEC, currentEC, maxScaleFactor);
    [actualNDisp, currentNDisp, nClim, nTicks, nLabels] = localLogDisplayData(actualN, currentN, maxScaleFactor);
    clayClim = localAlignedLinearClim(actualClay, currentClay, maxScaleFactor);

    nexttile; localPlotActualSlice(xDisplay, zDisplay, actualECDisp, ecClim, 'EC actual (analytic slice)', ecTicks, ecLabels);
    nexttile; localPlotCurrentSlice(mesh.x, mesh.z, currentECDisp, ecClim, 'EC current', ecTicks, ecLabels);
    nexttile; localPlotResidualSlice(mesh.x, mesh.z, residualECDisp, residualECClim, 'EC residual ((current-actual)/actual)', residualECTicks, residualECLabels);
    nexttile; plotVoltageFit(dObs, dPred);

    nexttile; localPlotActualSlice(xDisplay, zDisplay, actualNDisp, nClim, 'N actual (analytic slice)', nTicks, nLabels);
    nexttile; localPlotCurrentSlice(mesh.x, mesh.z, currentNDisp, nClim, 'N current', nTicks, nLabels);
    nexttile; localPlotResidualSlice(mesh.x, mesh.z, residualNDisp, residualNClim, 'N residual ((current-actual)/actual)', residualNTicks, residualNLabels);
    nexttile; plotMisfitConvergence(iterationLog, targetMisfit);

    nexttile; localPlotActualSlice(xDisplay, zDisplay, actualClay, clayClim, 'Clay actual (analytic slice)', [], []);
    nexttile; localPlotCurrentSlice(mesh.x, mesh.z, currentClay, clayClim, 'Clay current', [], []);
    nexttile; localPlotResidualSlice(mesh.x, mesh.z, residualClayDisp, residualClayClim, 'Clay residual ((current-actual)/actual)', residualClayTicks, residualClayLabels);
    nexttile; plotSummaryPanel(summary);
end

function localPlotActualSlice(xDisplay, zDisplay, values, climVals, ttl, cbTicks, cbLabels)
    imagesc(xDisplay, zDisplay, values);
    set(gca, 'YDir', 'reverse');
    axis tight;
    clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks)
        set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels);
    end
    colormap(gca, jet);
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function localPlotCurrentSlice(xEdges, zEdges, cellValues, climVals, ttl, cbTicks, cbLabels)
    plot_rectilinear_cells(xEdges, zEdges, cellValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks)
        set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels);
    end
    colormap(gca, jet);
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function localPlotResidualSlice(xEdges, zEdges, cellValues, climVals, ttl, cbTicks, cbLabels)
    plot_rectilinear_cells(xEdges, zEdges, cellValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    clim(climVals);
    cb = colorbar;
    if ~isempty(cbTicks)
        set(cb, 'Ticks', cbTicks, 'TickLabels', cbLabels);
    end
    colormap(gca, jet);
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function vals2 = localExtractXZSlice(mesh, values, yTarget)
    vals3 = reshape(values, mesh.gridSize - 1);
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - yTarget));
    vals2 = squeeze(vals3(:, iy, :))';
end

function sampled = localSampleActualToCoarse(mesh, actualGrid, xDisplay, zDisplay, yTarget)
    xCenters = 0.5 * (mesh.x(1:end-1) + mesh.x(2:end));
    zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - yTarget)); %#ok<ASGLU>
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
    if isempty(actualFinite) && isempty(currentFinite)
        climVals = [0, 1];
        return;
    end
    if isempty(actualFinite)
        actualFinite = currentFinite;
    end
    if isempty(currentFinite)
        currentFinite = actualFinite;
    end

    actualMin = min(actualFinite);
    actualMax = max(actualFinite);
    actualSpan = actualMax - actualMin;
    if ~(isfinite(actualSpan) && actualSpan > 0)
        actualSpan = max(abs(actualFinite));
    end
    if ~(isfinite(actualSpan) && actualSpan > 0)
        actualSpan = 1;
    end

    combinedMin = min([actualFinite; currentFinite]);
    combinedMax = max([actualFinite; currentFinite]);
    allowedSpan = maxScaleFactor * actualSpan;

    lower = combinedMin;
    upper = combinedMax;
    if (upper - lower) > allowedSpan
        lower = actualMin;
        upper = lower + allowedSpan;
        if combinedMax > upper
            upper = combinedMax;
            lower = upper - allowedSpan;
        end
    end

    if lower == upper
        pad = max(abs(lower), 1) * 0.5;
        climVals = [lower - pad, upper + pad];
    else
        climVals = [lower, upper];
    end
end

function plotVoltageFit(dObs, dPred)
    hold on;
    plot(dObs, dPred, 'o', 'Color', [0.35, 0.55, 0.95], ...
        'MarkerFaceColor', [0.35, 0.55, 0.95], ...
        'DisplayName', 'measurements');
    finiteVals = [dObs(:); dPred(:)];
    finiteVals = finiteVals(isfinite(finiteVals));
    if isempty(finiteVals)
        lims = [-1, 1];
    else
        lims = [min(finiteVals), max(finiteVals)];
        if lims(1) == lims(2)
            lims = lims + [-0.5, 0.5];
        end
    end
    plot(lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1');
    xlim(lims);
    ylim(lims);
    axis square;
    grid on;
    xlabel('observed voltage (V)');
    ylabel('current voltage (V)');
    title('Voltage fit');
    legend('Location', 'best');
end

function plotMisfitConvergence(iterationLog, targetMisfit)
    it = [iterationLog.iter];
    misfit = [iterationLog.normalizedDataMisfit];
    semilogy(it, max(misfit, realmin), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    hold on;
    yline(max(targetMisfit, realmin), 'r--', 'LineWidth', 1.3, 'DisplayName', 'target');
    grid on;
    set(gca, 'YScale', 'log');
    xlabel('iteration');
    ylabel('normalized misfit (log scale)');
    title('Residual convergence');
    legend({'current', 'target'}, 'Location', 'best');
end

function plotSummaryPanel(summary)
    axis off;
    x0 = 0.05;
    y = 0.90;
    dy = 0.11;
    text(x0, y, sprintf('Exit: %s', summary.exitReason), 'Units', 'normalized', 'FontWeight', 'bold');
    y = y - dy;
    text(x0, y, sprintf('Iterations: %d', summary.iterationsCompleted), 'Units', 'normalized');
    y = y - dy;
    text(x0, y, sprintf('Obj before: %.3g', summary.objectiveBefore), 'Units', 'normalized');
    y = y - dy;
    text(x0, y, sprintf('Obj after: %.3g', summary.objectiveAfter), 'Units', 'normalized');
    y = y - dy;
    text(x0, y, sprintf('Misfit after: %.3g', summary.normalizedDataMisfitAfter), 'Units', 'normalized');
    y = y - dy;
    if isfield(summary, 'nPriorRowsUsed') && summary.nPriorRowsUsed > 0
        text(x0, y, sprintf('Prior rows: %d', summary.nPriorRowsUsed), 'Units', 'normalized');
        y = y - dy;
    end
    if isfield(summary, 'priorPenaltyAfter')
        text(x0, y, sprintf('Prior penalty: %.3g', summary.priorPenaltyAfter), 'Units', 'normalized');
        y = y - dy;
    end
    if isfield(summary, 'priorEcPenaltyAfter') && summary.priorEcPenaltyAfter > 0
        text(x0, y, sprintf('Prior EC penalty: %.3g', summary.priorEcPenaltyAfter), 'Units', 'normalized');
        y = y - dy;
    end
    text(x0, y, sprintf('Step length: %.3g', summary.lastStepLength), 'Units', 'normalized');
    title('Run summary');
end

function [actualDisp, currentDisp, climVals, ticks, labels] = localLogDisplayData(actualVals, currentVals, maxScaleFactor)
    combined = [actualVals(:); currentVals(:)];
    positiveVals = combined(isfinite(combined) & combined > 0);
    if isempty(positiveVals)
        actualDisp = actualVals;
        currentDisp = currentVals;
        climVals = localAlignedLinearClim(actualVals, currentVals, maxScaleFactor);
        ticks = [];
        labels = {};
        return;
    end

    displayFloor = max(min(positiveVals), 1e-6);
    actualDisp = log10(max(actualVals, displayFloor));
    currentDisp = log10(max(currentVals, displayFloor));

    actualPositive = actualVals(isfinite(actualVals) & actualVals > 0);
    currentPositive = currentVals(isfinite(currentVals) & currentVals > 0);
    if isempty(actualPositive)
        actualPositive = positiveVals;
    end
    if isempty(currentPositive)
        currentPositive = actualPositive;
    end

    actualMin = max(min(actualPositive), displayFloor);
    actualMax = max(actualPositive);
    actualLogMin = log10(actualMin);
    actualLogMax = log10(actualMax);
    actualLogSpan = actualLogMax - actualLogMin;
    if ~(isfinite(actualLogSpan) && actualLogSpan > 0)
        actualLogSpan = 1;
    end

    combinedLogMin = min([log10(max(currentPositive, displayFloor)); actualLogMin]);
    combinedLogMax = max([log10(max(currentPositive, displayFloor)); actualLogMax]);
    allowedLogSpan = maxScaleFactor * actualLogSpan;

    lower = combinedLogMin;
    upper = combinedLogMax;
    if (upper - lower) > allowedLogSpan
        lower = actualLogMin;
        upper = lower + allowedLogSpan;
        if combinedLogMax > upper
            upper = combinedLogMax;
            lower = upper - allowedLogSpan;
        end
        lower = max(lower, log10(displayFloor));
    end

    if lower == upper
        climVals = [lower - 0.5, upper + 0.5];
    else
        climVals = [lower, upper];
    end

    visibleMin = 10 ^ climVals(1);
    visibleMax = 10 ^ climVals(2);
    if visibleMax <= visibleMin
        visibleMax = visibleMin * 10;
    end
    if numel(positiveVals) > 1
        tickVals = logspace(log10(visibleMin), log10(visibleMax), 4)';
    else
        tickVals = max(positiveVals(1), displayFloor);
    end
    tickVals = unique(tickVals);
    ticks = log10(tickVals);
    labels = arrayfun(@(x) sprintf('%.3g', x), tickVals, 'UniformOutput', false);
end

function [dispVals, climVals, ticks, labels] = localSignedLogResidualData(residualVals)
    finiteVals = residualVals(isfinite(residualVals));
    if isempty(finiteVals)
        dispVals = residualVals;
        climVals = [-1, 1];
        ticks = [];
        labels = {};
        return;
    end

    absVals = abs(finiteVals);
    absVals = absVals(absVals > 0);
    if isempty(absVals)
        dispVals = zeros(size(residualVals));
        climVals = [-1, 1];
        ticks = [-1, 0, 1];
        labels = {'-1', '0', '1'};
        return;
    end

    scale = max(prctile(absVals, 10), 1e-4);
    dispVals = sign(residualVals) .* log10(1 + abs(residualVals) ./ scale);

    maxMag = prctile(absVals, 98);
    if ~(isfinite(maxMag) && maxMag > 0)
        maxMag = max(absVals);
    end
    maxDisp = log10(1 + maxMag ./ scale);
    maxDisp = max(maxDisp, 1);
    climVals = [-maxDisp, maxDisp];

    magTicks = unique([scale; sqrt(scale * maxMag); maxMag]);
    posTicks = log10(1 + magTicks ./ scale);
    ticks = [-flipud(posTicks(:)); 0; posTicks(:)];
    labels = [arrayfun(@(x) sprintf('-%#.3g', x), flipud(magTicks(:)), 'UniformOutput', false); {'0'}; ...
              arrayfun(@(x) sprintf('%#.3g', x), magTicks(:), 'UniformOutput', false)];
end
