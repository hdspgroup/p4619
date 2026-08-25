function diagnose_direct_ec_random_forest(dataFile, modelFile)
% Diagnostic plots for the direct EC Random Forest model.
%
% For the expanded combined dataset:
%   diagnose_direct_ec_random_forest('data_comp.mat', 'direct_ec_rf_model.mat')

    if nargin < 1 || strlength(string(dataFile)) == 0
        dataFile = 'data_comp.mat';
    end
    if nargin < 2 || strlength(string(modelFile)) == 0
        modelFile = 'direct_ec_rf_model.mat';
    end

    paths = codex_paths('direct_ec_random_forest_data_comp');
    outDir = fullfile(paths.outputsDir, 'direct_ec_random_forest_diagnostics');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    dataPath = dataFile;
    if ~isfile(dataPath)
        dataPath = fullfile(paths.dataDir, dataFile);
    end

    modelPath = modelFile;
    if ~isfile(modelPath)
        modelPath = fullfile(paths.outputsDir, 'direct_ec_random_forest_data_comp', modelFile);
    end

    s = load(dataPath);
    m = load(modelPath);
    model = m.model;

    rawNames = string(s.varNames);
    names = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
    T = array2table(s.data, 'VariableNames', cellstr(names));

    targetName = model.metrics.TargetName;
    predictorNames = model.predictorNames;
    if targetName == "EC"
        T.EC = T.EC .* 0.1;
    end

    D = T(:, [predictorNames, targetName]);
    keep = all(~ismissing(D), 2);
    D = D(keep, :);

    Xraw = D{:, predictorNames};
    for v = model.skewedPredictors
        if ismember(v, predictorNames)
            D.(v) = log1p(D.(v));
        end
    end
    Xtrans = D{:, predictorNames};
    X = (Xtrans - model.mu) ./ model.sigma;
    y = D.(targetName);

    yhat = exp(predict(model.rf, X)) - 1;
    resid = yhat - y;
    logResid = log10(max(yhat, realmin)) - log10(max(y, realmin));

    qEdges = quantile(y, [0 1/3 2/3 1]);
    qEdges = unique(qEdges);
    if numel(qEdges) < 4
        qEdges = quantile(y, [0 0.25 0.5 0.75 1]);
        qEdges = unique(qEdges);
    end
    binIdx = discretize(y, qEdges, 'IncludedEdge', 'right');
    binLabels = strings(max(binIdx), 1);
    binRMSELog10 = nan(max(binIdx), 1);
    binMAE = nan(max(binIdx), 1);
    binMeanLogResid = nan(max(binIdx), 1);
    binMAELog10 = nan(max(binIdx), 1);
    for b = 1:max(binIdx)
        idx = binIdx == b;
        if ~any(idx)
            continue;
        end
        binLabels(b) = sprintf('[%.3g, %.3g]', qEdges(b), qEdges(b + 1));
        binRMSELog10(b) = sqrt(mean((logResid(idx)).^2));
        binMAE(b) = mean(abs(yhat(idx) - y(idx)));
        binMeanLogResid(b) = mean(logResid(idx));
        binMAELog10(b) = mean(abs(logResid(idx)));
    end
    binTbl = table(binLabels, binRMSELog10, binMAELog10, binMAE, binMeanLogResid, ...
        'VariableNames', {'ECBin','RMSE_log10','MAE_log10','MAE_Sm','MeanLogResidual'});
    writetable(binTbl, fullfile(outDir, 'ec_error_by_quantile.csv'));

    f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
    scatter(y, yhat, 10, 'filled', 'MarkerFaceAlpha', 0.2);
    hold on;
    lo = min([y; yhat]);
    hi = max([y; yhat]);
    plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Observed EC (S/m)');
    ylabel('Predicted EC (S/m)');
    title('Random Forest: observed vs predicted');
    grid on;
    exportgraphics(f1, fullfile(outDir, 'rf_observed_vs_predicted.png'), 'Resolution', 200);
    close(f1);

    f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
    scatter(y, logResid, 10, 'filled', 'MarkerFaceAlpha', 0.2);
    set(gca, 'XScale', 'log');
    yline(0, 'k--', 'LineWidth', 1.2);
    xlabel('Observed EC (S/m)');
    ylabel('log10(pred) - log10(obs)');
    title('Log residual vs observed EC');
    grid on;
    exportgraphics(f2, fullfile(outDir, 'rf_log_residual_vs_observed_ec.png'), 'Resolution', 200);
    close(f2);

    plotVars = intersect(["pH_CaCl2","pH_H2O","K","P","N","Clay","Silt","Sand","CaCO3","OC"], predictorNames, 'stable');
    nVars = numel(plotVars);
    nCols = 2;
    nRows = ceil(nVars / nCols);
    f3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 350*nRows]);
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:nVars
        nexttile;
        x = D.(plotVars(i));
        if any(model.skewedPredictors == plotVars(i))
            x = exp(x) - 1;
        end
        scatter(x, logResid, 8, 'filled', 'MarkerFaceAlpha', 0.2);
        yline(0, 'k--');
        xlabel(plotVars(i));
        ylabel('log residual');
        title(sprintf('Residual vs %s', plotVars(i)));
        grid on;
    end
    exportgraphics(f3, fullfile(outDir, 'rf_log_residual_vs_predictors.png'), 'Resolution', 180);
    close(f3);

    f4 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 500]);
    yyaxis left;
    bar(categorical(binTbl.ECBin), binTbl.RMSE_log10);
    ylabel('RMSE (log10 space)');
    yyaxis right;
    plot(categorical(binTbl.ECBin), binTbl.MeanLogResidual, '-o', 'LineWidth', 1.5);
    ylabel('Mean log residual');
    title('Error by EC quantile');
    grid on;
    exportgraphics(f4, fullfile(outDir, 'rf_error_by_ec_quantile.png'), 'Resolution', 200);
    close(f4);

    fid = fopen(fullfile(outDir, 'rf_diagnostic_summary.txt'), 'w');
    fprintf(fid, 'Random Forest diagnostic summary\n');
    fprintf(fid, 'Data file: %s\n', dataPath);
    fprintf(fid, 'Model file: %s\n', modelPath);
    fprintf(fid, 'Rows used: %d\n', numel(y));
    fprintf(fid, 'Overall RMSE = %.6g S/m\n', sqrt(mean((yhat - y).^2)));
    fprintf(fid, 'Overall MAE = %.6g S/m\n', mean(abs(yhat - y)));
    fprintf(fid, 'Overall mean log residual = %.6g\n\n', mean(logResid));
    fprintf(fid, 'Error by EC quantile:\n');
    for i = 1:height(binTbl)
        fprintf(fid, '  %s | RMSE_log10=%.6g | MAE_log10=%.6g | MAE_Sm=%.6g | mean log residual=%.6g\n', ...
            binTbl.ECBin(i), binTbl.RMSE_log10(i), binTbl.MAE_log10(i), binTbl.MAE_Sm(i), binTbl.MeanLogResidual(i));
    end
    fclose(fid);

    disp('Random Forest diagnostics complete.');
end
