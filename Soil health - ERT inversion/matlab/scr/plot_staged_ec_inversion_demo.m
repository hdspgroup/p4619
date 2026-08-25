function fig = plot_staged_ec_inversion_demo(mesh, survey, sigmaBackground, sigmaFull, ...
    trueSlices, dObs, dPredBackground, dPredFull, iterationLogBackground, iterationLogFull, summary)
%PLOT_STAGED_EC_INVERSION_DEMO Plot staged EC inversion results.

    ySlice = trueSlices.xz.y;
    zSlice = trueSlices.xy.z;

    xyActual = trueSlices.xy.v;
    xzActual = trueSlices.xz.v;
    xyBg = localExtractXYSlice(mesh, sigmaBackground, zSlice);
    xzBg = localExtractXZSlice(mesh, sigmaBackground, ySlice);
    xyFull = localExtractXYSlice(mesh, sigmaFull, zSlice);
    xzFull = localExtractXZSlice(mesh, sigmaFull, ySlice);

    ecRange = localAlignedRange([xyActual(:); xzActual(:)], [xyBg(:); xzBg(:)], [xyFull(:); xzFull(:)]);

    actualCoarseXY = localSampleActualToCoarse(mesh, xyActual, trueSlices.xy.x(1, :), trueSlices.xy.y(:, 1), 'xy');
    actualCoarseXZ = localSampleActualToCoarse(mesh, xzActual, trueSlices.xz.x(1, :), trueSlices.xz.z(:, 1), 'xz');
    xyResidual = localRelativeResidual(xyFull, actualCoarseXY);
    xzResidual = localRelativeResidual(xzFull, actualCoarseXZ);

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.02, 0.06, 0.96, 0.82]);
    tl = tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Staged EC inversion demo');

    nexttile(1); localPlotTrueXY(trueSlices.xy, ecRange, 'EC actual XY', survey);
    nexttile(2); localPlotCellXY(mesh, xyBg, zSlice, ecRange, 'Depth-background EC XY');
    nexttile(3); localPlotCellXY(mesh, xyFull, zSlice, ecRange, 'Background + anomaly EC XY');
    nexttile(4); localPlotResidualXY(mesh, xyResidual, zSlice, 'Residual XY ((full-actual)/actual)');
    nexttile(5); localPlotVoltageComparison(dObs, dPredBackground, dPredFull);

    nexttile(6); localPlotTrueXZ(trueSlices.xz, ecRange, 'EC actual XZ');
    nexttile(7); localPlotCellXZ(mesh, xzBg, ySlice, ecRange, 'Depth-background EC XZ');
    nexttile(8); localPlotCellXZ(mesh, xzFull, ySlice, ecRange, 'Background + anomaly EC XZ');
    nexttile(9); localPlotResidualXZ(mesh, xzResidual, ySlice, 'Residual XZ ((full-actual)/actual)');
    nexttile(10); localPlotConvergence(iterationLogBackground, iterationLogFull, summary);
end

function localPlotTrueXY(slice, climVals, ttl, survey)
    if isfield(slice, 'xEdges') && isfield(slice, 'yEdges')
        plot_rectilinear_cells(slice.xEdges, slice.yEdges, slice.v);
    else
        surface(slice.x, slice.y, zeros(size(slice.v)), slice.v, 'EdgeColor', 'none', 'FaceColor', 'flat');
        view(2);
    end
    axis equal tight;
    colormap(gca, jet); clim(climVals); colorbar;
    hold on;
    if ~isempty(survey) && isfield(survey, 'electrodes')
        plot(survey.electrodes(:, 1), survey.electrodes(:, 2), 'ko', 'MarkerSize', 4.5, 'LineWidth', 1.0);
    end
    hold off;
    title(sprintf('%s (z ~= %.1f m)', ttl, slice.z));
    xlabel('x (m)'); ylabel('y (m)');
end

function localPlotTrueXZ(slice, climVals, ttl)
    if isfield(slice, 'xEdges') && isfield(slice, 'zEdges')
        plot_rectilinear_cells(slice.xEdges, slice.zEdges, slice.v);
    else
        surface(slice.x, slice.z, zeros(size(slice.v)), slice.v, 'EdgeColor', 'none', 'FaceColor', 'flat');
        view(2);
    end
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('%s (y ~= %.1f m)', ttl, slice.y));
    xlabel('x (m)'); ylabel('z (m)');
end

function localPlotCellXY(mesh, vals2, zSlice, climVals, ttl)
    plot_rectilinear_cells(mesh.x, mesh.y, vals2);
    axis equal tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('%s (z ~= %.1f m)', ttl, zSlice));
    xlabel('x (m)'); ylabel('y (m)');
end

function localPlotCellXZ(mesh, vals2, ySlice, climVals, ttl)
    plot_rectilinear_cells(mesh.x, mesh.z, vals2);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('%s (y ~= %.1f m)', ttl, ySlice));
    xlabel('x (m)'); ylabel('z (m)');
end

function localPlotResidualXY(mesh, residualVals, zSlice, ttl)
    plot_rectilinear_cells(mesh.x, mesh.y, residualVals);
    axis equal tight;
    cmax = localResidualRange(residualVals);
    colormap(gca, jet); clim([-cmax, cmax]);
    cb = colorbar; ylabel(cb, '(full - actual) / actual');
    title(sprintf('%s (z ~= %.1f m)', ttl, zSlice));
    xlabel('x (m)'); ylabel('y (m)');
end

function localPlotResidualXZ(mesh, residualVals, ySlice, ttl)
    plot_rectilinear_cells(mesh.x, mesh.z, residualVals);
    set(gca, 'YDir', 'reverse'); axis tight;
    cmax = localResidualRange(residualVals);
    colormap(gca, jet); clim([-cmax, cmax]);
    cb = colorbar; ylabel(cb, '(full - actual) / actual');
    title(sprintf('%s (y ~= %.1f m)', ttl, ySlice));
    xlabel('x (m)'); ylabel('z (m)');
end

function localPlotVoltageComparison(dObs, dPredBackground, dPredFull)
    hold on;
    scatter(dObs, dPredBackground, 22, 'r', 'filled', 'MarkerFaceAlpha', 0.45, 'DisplayName', 'background');
    scatter(dObs, dPredFull, 22, 'b', 'filled', 'MarkerFaceAlpha', 0.55, 'DisplayName', 'full');
    lims = [min([dObs(:); dPredBackground(:); dPredFull(:)]), max([dObs(:); dPredBackground(:); dPredFull(:)])];
    if ~all(isfinite(lims)) || diff(lims) <= 0
        lims = [-1, 1];
    end
    plot(lims, lims, 'k--', 'LineWidth', 1.2, 'DisplayName', '1:1');
    xlim(lims); ylim(lims); axis square; grid on;
    xlabel('observed voltage (V)'); ylabel('current voltage (V)');
    title('Voltage fit');
    legend('Location', 'best');
end

function localPlotConvergence(iterationLogBackground, iterationLogFull, summary)
    hold on;
    if ~isempty(iterationLogBackground)
        itBg = 1:numel(iterationLogBackground);
        misBg = max([iterationLogBackground.normalizedDataMisfit], realmin);
        semilogy(itBg, misBg, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'background misfit');
    end
    if ~isempty(iterationLogFull)
        offset = numel(iterationLogBackground);
        itFull = offset + (1:numel(iterationLogFull));
        misFull = max([iterationLogFull.normalizedDataMisfit], realmin);
        semilogy(itFull, misFull, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'full misfit');
    end
    grid on; set(gca, 'YScale', 'log');
    xlabel('iteration'); ylabel('normalized misfit');
    title(sprintf('Convergence (%s)', summary.exitReason));
    legend('Location', 'best');
end

function vals2 = localExtractXYSlice(mesh, values, zSlice)
    vals3 = reshape(values, mesh.gridSize - 1);
    zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
    [~, iz] = min(abs(zCenters - zSlice));
    vals2 = squeeze(vals3(:, :, iz))';
end

function vals2 = localExtractXZSlice(mesh, values, ySlice)
    vals3 = reshape(values, mesh.gridSize - 1);
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - ySlice));
    vals2 = squeeze(vals3(:, iy, :))';
end

function sampled = localSampleActualToCoarse(mesh, actualGrid, axis1, axis2, mode)
    xCenters = 0.5 * (mesh.x(1:end-1) + mesh.x(2:end));
    if strcmp(mode, 'xy')
        yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
        [Xc, Yc] = meshgrid(xCenters, yCenters);
        sampled = interp2(axis1, axis2, actualGrid, Xc, Yc, 'linear');
    else
        zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));
        [Xc, Zc] = meshgrid(xCenters, zCenters);
        sampled = interp2(axis1, axis2, actualGrid, Xc, Zc, 'linear');
    end
end

function residual = localRelativeResidual(currentVals, actualVals)
    residual = (currentVals - actualVals) ./ max(abs(actualVals), 1e-6);
end

function vr = localAlignedRange(actualVals, bgVals, fullVals)
    vals = [actualVals(:); bgVals(:); fullVals(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vr = [0, 1];
        return;
    end
    vmin = min(vals); vmax = max(vals);
    if vmax <= vmin
        pad = max(abs(vmin), 1) * 0.05;
        vr = [vmin - pad, vmin + pad];
    else
        vr = [vmin, vmax];
    end
end

function cmax = localResidualRange(vals)
    finiteVals = abs(vals(isfinite(vals)));
    if isempty(finiteVals)
        cmax = 1;
        return;
    end
    cmax = max(prctile(finiteVals, 98), max(finiteVals) * 0.15);
    if ~(isfinite(cmax) && cmax > 0)
        cmax = 1;
    end
end
