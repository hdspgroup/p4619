function fig = plot_ec_only_inversion_demo(mesh, sigmaTrue, sigmaStart, sigmaCurrent, dObs, dPred, iterationLog, targetMisfit, summary, varargin)
%PLOT_EC_ONLY_INVERSION_DEMO Summary figure for EC-only inversion.

    p = inputParser;
    addParameter(p, 'ForwardMesh', []);
    addParameter(p, 'SigmaTrueForward', []);
    addParameter(p, 'Survey', []);
    addParameter(p, 'TolObjective', []);
    parse(p, varargin{:});
    forwardMesh = p.Results.ForwardMesh;
    sigmaTrueForward = p.Results.SigmaTrueForward;
    survey = p.Results.Survey;
    tolObjective = p.Results.TolObjective;

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.04, 0.08, 0.90, 0.78]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'EC-only inversion demo');

    yMid = round((mesh.gridSize(2) - 1) / 2);

    if ~isempty(forwardMesh) && ~isempty(sigmaTrueForward)
        ecRange = [ ...
            min([sigmaTrueForward(:); sigmaCurrent(:)]), ...
            max([sigmaTrueForward(:); sigmaCurrent(:)])];
    else
        ecRange = [ ...
            min([sigmaTrue(:); sigmaStart(:); sigmaCurrent(:)]), ...
            max([sigmaTrue(:); sigmaStart(:); sigmaCurrent(:)])];
    end
    if ~all(isfinite(ecRange)) || diff(ecRange) <= 0
        ecRange = [0, 1];
    end

    if ~isempty(forwardMesh) && ~isempty(sigmaTrueForward)
        yMidForward = round((forwardMesh.gridSize(2) - 1) / 2);
        nexttile; plotField(forwardMesh, sigmaTrueForward, yMidForward, 'EC actual fine (S/m)', ecRange);
        nexttile; plotField(mesh, sigmaCurrent, yMid, 'EC current coarse (S/m)', ecRange);
        nexttile; plotField(mesh, sigmaStart, yMid, 'EC start coarse (S/m)', ecRange);
    else
        nexttile; plotField(mesh, sigmaTrue, yMid, 'EC actual (S/m)', ecRange);
        nexttile; plotField(mesh, sigmaStart, yMid, 'EC start (S/m)', ecRange);
        nexttile; plotField(mesh, sigmaCurrent, yMid, 'EC current (S/m)', ecRange);
    end

    nexttile; plotCenteredResidualField(mesh, ...
        log10(max(sigmaCurrent, realmin)) - log10(max(sigmaTrue, realmin)), ...
        yMid, 'EC residual coarse (log_{10})');
    nexttile; plotVoltageFit(dObs, dPred, survey);
    nexttile; plotMisfitAndObjective(iterationLog, targetMisfit, tolObjective, summary);
end

function plotField(mesh, values, yMid, ttl, climVals)
    vals3 = reshape(values, mesh.gridSize - 1);
    cellVals = squeeze(vals3(:, yMid, :))';
    plot_rectilinear_cells(mesh.x, mesh.z, cellVals);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    colormap(gca, jet);
    if nargin >= 5 && ~isempty(climVals)
        clim(climVals);
    end
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotCenteredResidualField(mesh, values, yMid, ttl)
    vals3 = reshape(values, mesh.gridSize - 1);
    sliceVals = squeeze(vals3(:, yMid, :))';
    plot_rectilinear_cells(mesh.x, mesh.z, sliceVals);
    set(gca, 'YDir', 'reverse');
    axis tight;
    colormap(gca, jet);
    cmax = max(abs(sliceVals(:)));
    if ~isfinite(cmax) || cmax == 0
        cmax = 1;
    end
    clim([-cmax, cmax]);
    cb = colorbar;
    ylabel(cb, 'current - actual');
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotVoltageFit(dObs, dPred, survey)
    if nargin >= 3 && ~isempty(survey) && isfield(survey, 'quads') && isfield(survey, 'electrodes')
        xVals = localPotentialPairMidpoints(survey);
        [xVals, order] = sort(xVals);
        dObs = dObs(order);
        dPred = dPred(order);

        hold on;
        plot(xVals, dObs, 'k-o', 'DisplayName', 'observed');
        plot(xVals, dPred, 'r-^', 'DisplayName', 'current');
        xlabel('M-N midpoint x (m)');
    else
        hold on;
        plot(dObs, 'k-o', 'DisplayName', 'observed');
        plot(dPred, 'r-^', 'DisplayName', 'current');
        xlabel('measurement');
    end
    grid on;
    ylabel('voltage (V)');
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
    grid on;
    xlabel('iteration');
    ylabel('log scale');
    title(sprintf('Convergence (%s)', summary.exitReason));
end

function xMid = localPotentialPairMidpoints(survey)
%LOCALPOTENTIALPAIRMIDPOINTS Midpoint x-coordinate of each M-N pair.

    quads = survey.quads;
    elec = survey.electrodes;
    mIds = quads(:, 3);
    nIds = quads(:, 4);
    xMid = 0.5 * (elec(mIds, 1) + elec(nIds, 1));
end
