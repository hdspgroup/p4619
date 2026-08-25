function fig = plot_ert_forward_results(mesh, survey, result, options)
%PLOT_ERT_FORWARD_RESULTS Plot surface and cross-section potential fields.

    arguments
        mesh struct
        survey struct
        result struct
        options.ReferenceResult = []
        options.PropertyElements (:,1) double = []
        options.Anomaly struct = struct()
        options.SignedLogMinAbs (1,1) double = 1e-3
    end

    phi = result.phiBySource(:, 1);
    phi3 = reshape(phi, mesh.gridSize - 1);

    if ~isempty(options.ReferenceResult)
        phiRef = options.ReferenceResult.phiBySource(:, 1);
        phiRef3 = reshape(phiRef, mesh.gridSize - 1);
        rel = (phi3 - phiRef3) ./ max(abs(phiRef3), realmin);
    else
        rel = [];
    end

    yMid = round((numel(mesh.y) - 1) / 2);
    zSurf = 1;

    fig = figure('Color', 'w', 'units', 'normalized',...
        'Position', [0.1, 0.05, 0.8, 0.8]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot_surface_panel(mesh, survey, phi3(:, :, zSurf), ...
        options.Anomaly, options.SignedLogMinAbs);
    title('surface \phi');

    nexttile;
    plot_xz_panel(mesh, survey, squeeze(phi3(:, yMid, :)), ...
        options.Anomaly, options.SignedLogMinAbs);
    title('\phi at mid Y');

    nexttile;
    if isempty(rel)
        plot_property_surface(mesh, options.PropertyElements, options.Anomaly);
        title('surface resistivity');
        cb = colorbar;
        ylabel(cb, '\rho (\Omega m)');
    else
        plot_surface_panel(mesh, survey, rel(:, :, zSurf), ...
            options.Anomaly, options.SignedLogMinAbs);
        title('relative potential difference');
    end

    nexttile;
    if isempty(rel)
        plot_property_xz(mesh, options.PropertyElements, options.Anomaly, yMid);
        title('resistivity at mid Y');
        cb = colorbar;
        ylabel(cb, '\rho (\Omega m)');
    else
        plot_xz_panel(mesh, survey, squeeze(rel(:, yMid, :)), ...
            options.Anomaly, options.SignedLogMinAbs);
        title('relative difference at mid Y');
    end
end

function plot_surface_panel(mesh, survey, values, anomaly, minAbs)
    [plotValues, ticks, labels] = signed_log_color_data(values', minAbs);
    plot_rectilinear_cells(mesh.x, mesh.y, plotValues);
    axis equal tight;
    colormap jet;
    clim([0, 1]);
    cb = colorbar;
    set(cb, 'Ticks', ticks, 'TickLabels', labels);
    ylabel(cb, 'signed log_{10} magnitude');
    hold on;
    plot_electrodes(survey);
    plot_anomaly_surface(anomaly);
    xlabel('x (m)');
    ylabel('y (m)');
end

function plot_xz_panel(mesh, survey, values, anomaly, minAbs)
    [plotValues, ticks, labels] = signed_log_color_data(values', minAbs);
    plot_rectilinear_cells(mesh.x, mesh.z, plotValues);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap jet;
    clim([0, 1]);
    cb = colorbar;
    set(cb, 'Ticks', ticks, 'TickLabels', labels);
    ylabel(cb, 'signed log_{10} magnitude');
    hold on;
    plot(survey.electrodes(:, 1), survey.electrodes(:, 3), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5);
    plot_anomaly_xz(anomaly);
    xlabel('x (m)');
    ylabel('z (m)');
end

function [plotValues, ticks, labels] = signed_log_color_data(values, minAbs)
    if nargin < 2 || isempty(minAbs) || minAbs <= 0
        minAbs = 1e-3;
    end

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

    tickExponents = logMin:ceil(max(1, (logMax-logMin)/4)):logMax;
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

function plot_property_surface(mesh, propertyElements, anomaly)
    if isempty(propertyElements)
        axis off;
        text(0.5, 0.5, 'No property field supplied', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        return;
    end
    values = element_values_to_grid(mesh, propertyElements);
    plot_rectilinear_cells(mesh.x, mesh.y, values(:, :, 1)');
    axis equal tight;
    colormap jet;
    hold on;
    plot_anomaly_surface(anomaly);
    xlabel('x (m)');
    ylabel('y (m)');
end

function plot_property_xz(mesh, propertyElements, anomaly, yMid)
    if isempty(propertyElements)
        axis off;
        text(0.5, 0.5, 'No property field supplied', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        return;
    end
    values = element_values_to_grid(mesh, propertyElements);
    yMidElement = min(yMid, size(values, 2));
    plot_rectilinear_cells(mesh.x, mesh.z, squeeze(values(:, yMidElement, :))');
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap jet;
    hold on;
    plot_anomaly_xz(anomaly);
    xlabel('x (m)');
    ylabel('z (m)');
end

function values = element_values_to_grid(mesh, propertyElements)
    n = mesh.gridSize - 1;
    values = reshape(propertyElements, n);
end

function plot_electrodes(survey)
    cols = {'r', 'b', 'g', 'g'};
    for i = 1:min(4, survey.nElectrodes)
        plot(survey.electrodes(i, 1), survey.electrodes(i, 2), '.', ...
            'Color', cols{i}, 'MarkerSize', 22);
    end
end

function plot_anomaly_surface(anomaly)
    if ~isfield(anomaly, 'center') || ~isfield(anomaly, 'radius')
        return;
    end
    t = linspace(0, 2*pi, 160);
    x = anomaly.center(1) + anomaly.radius * cos(t);
    y = anomaly.center(2) + anomaly.radius * sin(t);
    plot(x, y, 'k--', 'LineWidth', 1);
end

function plot_anomaly_xz(anomaly)
    if ~isfield(anomaly, 'center') || ~isfield(anomaly, 'radius')
        return;
    end
    t = linspace(0, 2*pi, 160);
    x = anomaly.center(1) + anomaly.radius * cos(t);
    z = anomaly.center(3) + anomaly.radius * sin(t);
    plot(x, z, 'k--', 'LineWidth', 1);
end
