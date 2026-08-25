function fig = plot_property_fit_ratios(mesh, fieldsTrue, fieldsCurrent, sigmaTrue, sigmaCurrent)
%PLOT_PROPERTY_FIT_RATIOS Plot model/actual ratios with a log-scaled colorbar.

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'normalized', ...
        'Position', [0.05, 0.06, 0.90, 0.82]);
    tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Property fit ratios: model / actual');

    yMid = round((mesh.gridSize(2) - 1) / 2);

    nexttile; plotRatioField(mesh, sigmaCurrent, sigmaTrue, yMid, 'EC ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.N, fieldsTrue.N, yMid, 'N ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.P, fieldsTrue.P, yMid, 'P ratio');

    nexttile; plotRatioField(mesh, fieldsCurrent.K, fieldsTrue.K, yMid, 'K ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.OC, fieldsTrue.OC, yMid, 'OC ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.Moisture, fieldsTrue.Moisture, yMid, 'Moisture ratio');

    nexttile; plotRatioField(mesh, fieldsCurrent.Clay, fieldsTrue.Clay, yMid, 'Clay ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.Sand, fieldsTrue.Sand, yMid, 'Sand ratio');
    nexttile; plotRatioField(mesh, fieldsCurrent.Silt, fieldsTrue.Silt, yMid, 'Silt ratio');
end

function plotRatioField(mesh, modelValues, trueValues, yMid, ttl)
    ratio = safe_ratio(modelValues, trueValues);
    vals3 = reshape(ratio, mesh.gridSize - 1);
    ratioSlice = squeeze(vals3(:, yMid, :))';

    [plotValues, ticks, labels] = ratio_color_data(ratioSlice);
    imagesc(mesh.x(1:end-1), mesh.z(1:end-1), plotValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap jet;
    clim([0, 1]);
    cb = colorbar;
    set(cb, 'Ticks', ticks, 'TickLabels', labels);
    ylabel(cb, 'model / actual');
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function [plotValues, ticks, labels] = ratio_color_data(ratioValues)
    finiteVals = ratioValues(isfinite(ratioValues) & ratioValues > 0);
    plotValues = 0.5 * ones(size(ratioValues));

    if isempty(finiteVals)
        ticks = 0.5;
        labels = {'1'};
        return;
    end

    logRatio = log10(finiteVals);
    logAbsMax = max(abs(logRatio));
    logAbsMax = min(max(logAbsMax, log10(1.25)), 2);

    fullLogRatio = log10(max(ratioValues, realmin));
    scaled = 0.5 + 0.5 * (fullLogRatio / logAbsMax);
    plotValues(isfinite(scaled)) = max(0, min(1, scaled(isfinite(scaled))));

    tickExponents = unique([-logAbsMax, -logAbsMax/2, 0, logAbsMax/2, logAbsMax]);
    ticks = 0.5 + 0.5 * (tickExponents / logAbsMax);
    ratioTicks = 10 .^ tickExponents;
    labels = arrayfun(@(x) sprintf('%.2g', x), ratioTicks, 'UniformOutput', false);
end

function ratio = safe_ratio(modelValues, trueValues)
    tiny = 1e-6;
    maxRatio = 1e2;
    minRatio = 1 / maxRatio;

    ratio = ones(size(modelValues));
    trueSmall = abs(trueValues) < tiny;
    modelSmall = abs(modelValues) < tiny;

    normalMask = ~trueSmall;
    ratio(normalMask) = max(modelValues(normalMask), tiny) ./ max(trueValues(normalMask), tiny);

    overMask = trueSmall & ~modelSmall;
    ratio(overMask) = maxRatio;

    underMask = ~trueSmall & modelSmall;
    ratio(underMask) = minRatio;

    bothSmall = trueSmall & modelSmall;
    ratio(bothSmall) = 1;

    ratio = min(max(ratio, minRatio), maxRatio);
end
