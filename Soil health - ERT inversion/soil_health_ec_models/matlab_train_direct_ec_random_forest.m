function model = train_direct_ec_random_forest(dataFile, outputModelFile)
% Train a direct Random Forest model for EC_sat from current soil variables.
%
% Example:
%   model = train_direct_ec_random_forest('data.mat', 'direct_ec_rf_model.mat');

    if nargin < 1 || strlength(string(dataFile)) == 0
        dataFile = 'data.mat';
    end
    paths = codex_paths('direct_ec_random_forest');
    outDir = paths.analysisDir;
    if ~isfile(dataFile)
        dataFile = fullfile(paths.dataDir, dataFile);
    end
    if nargin < 2 || strlength(string(outputModelFile)) == 0
        outputModelFile = fullfile(outDir, 'direct_ec_rf_model.mat');
    elseif ~contains(char(outputModelFile), filesep)
        outputModelFile = fullfile(outDir, char(outputModelFile));
    end

    s = load(dataFile);
    rawNames = string(s.varNames);
    names = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
    T = array2table(s.data, 'VariableNames', cellstr(names));

    if ismember("EC_sat", string(T.Properties.VariableNames))
        predictorNames = ["pH_CaCl2","pH_H2O","OC","P","N","K","Coarse","Clay","Sand","Silt"];
        targetName = "EC_sat";
        skewedPredictors = ["OC","P","N","K","Coarse","Clay"];
        analysisName = 'direct_ec_random_forest';
        targetUnits = 'S/m';
    elseif ismember("EC", string(T.Properties.VariableNames))
        predictorNames = ["Coarse","Clay","Silt","Sand","pH_CaCl2","pH_H2O","OC","CaCO3","N","P","K"];
        targetName = "EC";
        skewedPredictors = ["OC","CaCO3","N","P","K","Coarse","Clay"];
        analysisName = 'direct_ec_random_forest_data_comp';
        targetUnits = 'dS/m input, converted to S/m';
        T.EC = T.EC .* 0.1;
    else
        error('No supported EC target found. Expected EC_sat or EC.');
    end

    outDir = fullfile(paths.outputsDir, analysisName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    D = T(:, [predictorNames, targetName]);
    keep = all(~ismissing(D), 2);
    D = D(keep, :);

    [X, y, mu, sigma] = localPrepareData(D, predictorNames, skewedPredictors, targetName);

    rf = TreeBagger(400, X, log1p(y), ...
        'Method', 'regression', ...
        'MinLeafSize', 4, ...
        'NumPredictorsToSample', max(1, round(sqrt(numel(predictorNames)))), ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on');

    yhat = exp(predict(rf, X)) - 1;
    [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat);

    metrics = struct( ...
        'N', numel(y), ...
        'RMSE', rmse, ...
        'MAE', mae, ...
        'R2_log10', r2, ...
        'RSS_log10', rss, ...
        'TSS_log10', tss, ...
        'TargetName', targetName, ...
        'TargetUnits', targetUnits);

    importanceTbl = table(predictorNames', rf.OOBPermutedPredictorDeltaError(:), ...
        'VariableNames', {'Predictor','OOBImportance'});
    importanceTbl = sortrows(importanceTbl, 'OOBImportance', 'descend');

    model = struct();
    model.rf = rf;
    model.predictorNames = predictorNames;
    model.skewedPredictors = skewedPredictors;
    model.mu = mu;
    model.sigma = sigma;
    model.metrics = metrics;
    model.importance = importanceTbl;

    save(outputModelFile, 'model');
    writetable(importanceTbl, fullfile(outDir, 'direct_ec_rf_importance.csv'));

    fid = fopen(fullfile(outDir, 'direct_ec_rf_summary.txt'), 'w');
    fprintf(fid, 'Direct EC Random Forest summary\n');
    fprintf(fid, 'Target: %s\n', metrics.TargetName);
    fprintf(fid, 'Units: %s\n', metrics.TargetUnits);
    fprintf(fid, 'Training rows: %d\n', metrics.N);
    fprintf(fid, 'RMSE = %.6g\n', metrics.RMSE);
    fprintf(fid, 'MAE = %.6g\n', metrics.MAE);
    fprintf(fid, 'RSS in log10 space = %.6g\n', metrics.RSS_log10);
    fprintf(fid, 'TSS in log10 space = %.6g\n', metrics.TSS_log10);
    fprintf(fid, 'R^2 in log10 space = %.4f\n', metrics.R2_log10);
    fclose(fid);

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 820 700]);
    obs = max(y, realmin);
    pred = max(yhat, realmin);
    scatter(obs, pred, 12, 'filled', 'MarkerFaceAlpha', 0.28);
    hold on;
    lo = min([obs; pred]);
    hi = max([obs; pred]);
    plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Actual EC_{sat}');
    ylabel('Predicted EC_{sat}');
    title('Direct EC Random Forest');
    txt = sprintf('R^2_{log10} = %.3f | RMSE = %.3g | MAE = %.3g', ...
        metrics.R2_log10, metrics.RMSE, metrics.MAE);
    text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 6);
    grid on;
    exportgraphics(f, fullfile(outDir, 'direct_ec_rf_predicted_vs_actual.png'), 'Resolution', 200);
    close(f);
end

function [X, y, mu, sigma] = localPrepareData(D, predictorNames, skewedPredictors, targetName)
    for v = skewedPredictors
        D.(v) = log1p(D.(v));
    end
    Xraw = D{:, predictorNames};
    y = D.(targetName);
    mu = mean(Xraw, 1);
    sigma = std(Xraw, 0, 1);
    sigma(sigma == 0) = 1;
    X = (Xraw - mu) ./ sigma;
end

function [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat)
    keep = isfinite(y) & isfinite(yhat);
    y = y(keep);
    yhat = yhat(keep);
    rmse = sqrt(mean((y - yhat).^2));
    mae = mean(abs(y - yhat));
    yLog = log10(max(y, realmin));
    yhatLog = log10(max(yhat, realmin));
    tss = sum((yLog - mean(yLog)).^2);
    rss = sum((yLog - yhatLog).^2);
    r2 = 1 - rss / tss;
end
