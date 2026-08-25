function fig = plot_ert_pseudosection(survey, dObs, dPred)
%PLOT_ERT_PSEUDOSECTION Plot observed, predicted, and residual data by array geometry.

    [xMid, zPseudo] = pseudosection_coordinates(survey);
    resid = dPred - dObs;

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.08, 0.10, 0.82, 0.70]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'ERT pseudosection');

    nexttile;
    plotPseudoPanel(xMid, zPseudo, dObs, 'Observed voltage (V)', false);

    nexttile;
    plotPseudoPanel(xMid, zPseudo, dPred, 'Current voltage (V)', false);

    nexttile;
    plotPseudoPanel(xMid, zPseudo, resid, 'Residual (current - observed)', true);
end

function [xMid, zPseudo] = pseudosection_coordinates(survey)
    quads = survey.quads;
    elec = survey.electrodes;
    n = size(quads, 1);
    xMid = zeros(n, 1);
    zPseudo = zeros(n, 1);

    for i = 1:n
        ids = quads(i, :);
        x = elec(ids, 1);
        xMid(i) = mean(x);
        span = max(x) - min(x);
        zPseudo(i) = span / 6;
    end
end

function plotPseudoPanel(xMid, zPseudo, values, ttl, centerZero)
    scatterSize = 180;
    if centerZero
        vmax = max(abs(values));
        if ~isfinite(vmax) || vmax == 0
            vmax = 1;
        end
        scatter(xMid, zPseudo, scatterSize, values, 's', 'filled');
        colormap(gca, jet);
        clim([-vmax, vmax]);
        cb = colorbar;
        ylabel(cb, ttl);
    else
        scatter(xMid, zPseudo, scatterSize, values, 's', 'filled');
        colormap(gca, jet);
        cb = colorbar;
        ylabel(cb, ttl);
    end
    set(gca, 'YDir', 'reverse');
    grid on;
    xlabel('x midpoint (m)');
    ylabel('pseudodepth (m)');
    title(ttl);
end
