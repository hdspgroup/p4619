function [fig, selection] = plot_selected_source_potential(forwardMesh, survey, resultTrue, resultCurrent, options)
%PLOT_SELECTED_SOURCE_POTENTIAL Plot actual/current potential field for one A-B pair.

    arguments
        forwardMesh struct
        survey struct
        resultTrue struct
        resultCurrent struct
        options.SourceIndex (1,1) double = NaN
        options.SignedLogMinAbs (1,1) double = 1e-3
    end

    sourcePairs = resultTrue.sourcePairs;
    elec = survey.electrodes;

    if isnan(options.SourceIndex)
        idx = localChooseSourcePair(sourcePairs, elec);
    else
        idx = options.SourceIndex;
    end

    pair = sourcePairs(idx, :);
    aX = elec(pair(1), 1);
    bX = elec(pair(2), 1);
    pairMid = 0.5 * (aX + bX);
    pairSpan = abs(bX - aX);

    phiTrue = resultTrue.phiBySource(:, idx);
    phiCurrent = resultCurrent.phiBySource(:, idx);
    phiResidual = phiCurrent - phiTrue;

    phiTrue3 = reshape(phiTrue, forwardMesh.gridSize - 1);
    phiCurrent3 = reshape(phiCurrent, forwardMesh.gridSize - 1);
    phiResidual3 = reshape(phiResidual, forwardMesh.gridSize - 1);
    yMid = round((numel(forwardMesh.y) - 1) / 2);

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05, 0.08, 0.88, 0.78]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Potential field for A-B pair [%d, %d]  |  midpoint = %.1f m, span = %.1f m', ...
        pair(1), pair(2), pairMid, pairSpan));

    nexttile;
    localPlotPotentialXZ(forwardMesh, survey, squeeze(phiTrue3(:, yMid, :)), pair, ...
        options.SignedLogMinAbs, 'Actual \phi at mid Y');

    nexttile;
    localPlotPotentialXZ(forwardMesh, survey, squeeze(phiCurrent3(:, yMid, :)), pair, ...
        options.SignedLogMinAbs, 'Current \phi at mid Y');

    nexttile;
    localPlotResidualXZ(forwardMesh, survey, squeeze(phiResidual3(:, yMid, :)), pair, ...
        'Residual \phi (current - actual)');

    selection = struct();
    selection.sourceIndex = idx;
    selection.pair = pair;
    selection.midpoint = pairMid;
    selection.span = pairSpan;
end

function idx = localChooseSourcePair(sourcePairs, electrodes)
%LOCALCHOOSESOURCEPAIR Prefer a central pair with a large A-B span.

    aX = electrodes(sourcePairs(:, 1), 1);
    bX = electrodes(sourcePairs(:, 2), 1);
    mid = 0.5 * (aX + bX);
    span = abs(bX - aX);

    midNorm = abs(mid) / max(max(abs(mid)), eps);
    spanNorm = span / max(max(span), eps);
    score = midNorm - 0.75 * spanNorm;
    [~, idx] = min(score);
end

function localPlotPotentialXZ(mesh, survey, values, pair, minAbs, ttl)
    cellValues = localSliceToCellValues(mesh, values);
    [plotValues, ticks, labels] = localSignedLogColorData(cellValues, minAbs);
    plot_rectilinear_cells(mesh.x, mesh.z, plotValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap(gca, jet);
    clim([0, 1]);
    cb = colorbar;
    set(cb, 'Ticks', ticks, 'TickLabels', labels);
    ylabel(cb, 'signed log_{10} magnitude');
    hold on;
    localPlotSelectedElectrodes(survey, pair);
    xlabel('x (m)');
    ylabel('z (m)');
    title(ttl);
end

function localPlotResidualXZ(mesh, survey, values, pair, ttl)
    cellValues = localSliceToCellValues(mesh, values);
    plot_rectilinear_cells(mesh.x, mesh.z, cellValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap(gca, jet);
    cmax = max(abs(values(:)));
    if ~isfinite(cmax) || cmax == 0
        cmax = 1;
    end
    clim([-cmax, cmax]);
    cb = colorbar;
    ylabel(cb, '\phi residual (V)');
    hold on;
    localPlotSelectedElectrodes(survey, pair);
    xlabel('x (m)');
    ylabel('z (m)');
    title(ttl);
end

function localPlotSelectedElectrodes(survey, pair)
    plot(survey.electrodes(:, 1), survey.electrodes(:, 3), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4);
    plot(survey.electrodes(pair(1), 1), survey.electrodes(pair(1), 3), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 7);
    plot(survey.electrodes(pair(2), 1), survey.electrodes(pair(2), 3), 'bo', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 7);
end

function [plotValues, ticks, labels] = localSignedLogColorData(values, minAbs)
    finiteVals = values(isfinite(values) & abs(values) >= minAbs);
    plotValues = 0.5 * ones(size(values));

    if isempty(finiteVals)
        ticks = 0.5;
        labels = {'0'};
        return;
    end

    logAbs = log10(abs(finiteVals));
    logMax = ceil(max(logAbs));
    logMin = log10(minAbs);
    if logMax == logMin
        logMin = logMin - 1;
    end

    fullLogAbs = log10(abs(values));
    mag = (fullLogAbs - logMin) ./ (logMax - logMin);
    mag = max(0, min(1, mag));

    pos = values > 0 & isfinite(values) & abs(values) >= minAbs;
    neg = values < 0 & isfinite(values) & abs(values) >= minAbs;
    plotValues(pos) = 0.5 + 0.5 * mag(pos);
    plotValues(neg) = 0.5 - 0.5 * mag(neg);

    tickExponents = logMin:ceil(max(1, (logMax - logMin) / 4)):logMax;
    tickExponents = unique([logMin, tickExponents, logMax]);
    sideExponents = tickExponents(tickExponents > logMin);
    if isempty(sideExponents)
        sideExponents = logMax;
    end
    magTicks = (sideExponents - logMin) ./ (logMax - logMin);

    posTicks = 0.5 + 0.5 * magTicks;
    negTicks = 0.5 - 0.5 * fliplr(magTicks);
    ticks = [negTicks, 0.5, posTicks];

    negLabels = strings(1, numel(sideExponents));
    posLabels = strings(1, numel(sideExponents));
    revExp = fliplr(sideExponents);
    for i = 1:numel(sideExponents)
        negLabels(i) = sprintf('-10^{%d}', revExp(i));
        posLabels(i) = sprintf('10^{%d}', sideExponents(i));
    end
    labels = cellstr([negLabels, "0", posLabels]);

    [ticks, ia] = unique(ticks, 'stable');
    labels = labels(ia);
end

function cellValues = localSliceToCellValues(mesh, values)
    if ndims(values) ~= 2
        error('values must be a 2D x-z slice.');
    end

    nx = numel(mesh.x) - 1;
    nz = numel(mesh.z) - 1;
    if isequal(size(values), [nx, nz])
        cellValues = values';
        return;
    end
    if isequal(size(values), [numel(mesh.x), numel(mesh.z)])
        cellValuesXZ = 0.25 * ( ...
            values(1:end-1, 1:end-1) + ...
            values(2:end,   1:end-1) + ...
            values(1:end-1, 2:end) + ...
            values(2:end,   2:end));
        cellValues = cellValuesXZ';
        return;
    end

    error('Slice size does not match either cell-centered or nodal x-z dimensions.');
end
