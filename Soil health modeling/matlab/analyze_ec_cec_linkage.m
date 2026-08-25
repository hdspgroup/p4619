function analyze_ec_cec_linkage()
% Compare ways to use LUCAS repeat-point EC and CEC data.
%
% EC is converted from dS/m to S/m.
%
% Outputs:
%   ec_cec_linkage_summary.txt
%   ec_cec_model_comparison.csv
%   stable_point_table.csv
%   stable_cec_ec_scatter.png

    close all force;
    clc;

    paths = codex_paths('ec_cec_linkage');
    outDir = paths.analysisDir;

    s = load(fullfile(paths.dataDir, 'data_comp.mat'));
    varNames = matlab.lang.makeValidName(string(s.varNames), 'ReplacementStyle', 'delete');
    rawTable = array2table(s.data, 'VariableNames', cellstr(varNames));
    stableRowTable = array2table(s.stableData, 'VariableNames', cellstr(varNames));

    rawTable.EC = rawTable.EC .* 0.1;
    stableRowTable.EC = stableRowTable.EC .* 0.1;

    predictorBase = ["Coarse","Clay","Silt","Sand","pH_CaCl2","pH_H2O","OC","CaCO3","N","P","K"];
    targetEC = "EC";
    targetCEC = "CEC";

    stablePointTable = localBuildStablePointTable(stableRowTable, s.stableflag);
    writetable(stablePointTable, fullfile(outDir, 'stable_point_table.csv'));

    rawEC = localComplete(rawTable, [predictorBase, targetEC]);
    rawCEC = localComplete(rawTable, [predictorBase, targetCEC]);
    stableBase = localComplete(stablePointTable, [predictorBase, targetEC]);
    stableWithCEC = localComplete(stablePointTable, [predictorBase, targetCEC, targetEC]);

    rawECMetrics = localCrossValidatedRF(rawEC, predictorBase, targetEC, false);
    rawCECMetrics = localCrossValidatedRF(rawCEC, predictorBase, targetCEC, false);
    stableBaseMetrics = localCrossValidatedRF(stableBase, predictorBase, targetEC, false);
    stableCECMetrics = localCrossValidatedRF(stableWithCEC, [predictorBase, targetCEC], targetEC, false);

    comparison = table( ...
        ["Raw EC rows, no CEC"; "Raw CEC rows, CEC target"; "Stable point EC, no CEC"; "Stable point EC, measured/linked CEC"], ...
        [height(rawEC); height(rawCEC); height(stableBase); height(stableWithCEC)], ...
        [rawECMetrics.R2_log10; rawCECMetrics.R2_log10; stableBaseMetrics.R2_log10; stableCECMetrics.R2_log10], ...
        [rawECMetrics.RMSE; rawCECMetrics.RMSE; stableBaseMetrics.RMSE; stableCECMetrics.RMSE], ...
        [rawECMetrics.MAE; rawCECMetrics.MAE; stableBaseMetrics.MAE; stableCECMetrics.MAE], ...
        'VariableNames', {'Model','N','R2_log10','RMSE','MAE'});
    writetable(comparison, fullfile(outDir, 'ec_cec_model_comparison.csv'));

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 650]);
    keep = isfinite(stablePointTable.CEC) & isfinite(stablePointTable.EC) & stablePointTable.EC > 0;
    scatter(stablePointTable.CEC(keep), stablePointTable.EC(keep), 18, 'filled', 'MarkerFaceAlpha', 0.35);
    set(gca, 'YScale', 'log');
    xlabel('Linked CEC');
    ylabel('EC (S/m)');
    title('Stable repeat points: CEC vs EC');
    grid on;
    exportgraphics(f, fullfile(outDir, 'stable_cec_ec_scatter.png'), 'Resolution', 200);
    close(f);

    fid = fopen(fullfile(outDir, 'ec_cec_linkage_summary.txt'), 'w');
    fprintf(fid, 'EC-CEC linkage analysis\n');
    fprintf(fid, 'Rows in original data: %d\n', height(rawTable));
    fprintf(fid, 'Rows with EC: %d\n', sum(isfinite(rawTable.EC)));
    fprintf(fid, 'Rows with CEC: %d\n', sum(isfinite(rawTable.CEC)));
    fprintf(fid, 'Rows with both EC and CEC: %d\n', sum(isfinite(rawTable.EC) & isfinite(rawTable.CEC)));
    fprintf(fid, 'Stable repeat-point groups: %d\n', height(stablePointTable));
    fprintf(fid, 'Stable groups with EC and CEC linked by POINTID: %d\n\n', height(stableWithCEC));
    fprintf(fid, 'Model comparison, EC in S/m, R^2 in log10 target space:\n');
    for i = 1:height(comparison)
        fprintf(fid, '  %s: N=%d, R^2=%.4f, RMSE=%.6g, MAE=%.6g\n', ...
            comparison.Model(i), comparison.N(i), comparison.R2_log10(i), comparison.RMSE(i), comparison.MAE(i));
    end
    fprintf(fid, '\nInterpretation guidance:\n');
    fprintf(fid, '  Use the raw EC model for the strongest direct predictive benchmark.\n');
    fprintf(fid, '  Use the stable point linked model to test whether CEC adds information at repeated stable sites.\n');
    fprintf(fid, '  Do not treat linked CEC as same-year measured CEC; it is a point-level proxy transferred across campaigns.\n');
    fclose(fid);

    disp('EC-CEC linkage analysis complete.');
end

function pointTable = localBuildStablePointTable(T, stableflag)
    varNames = string(T.Properties.VariableNames);
    rows = stableflag;
    nGroups = size(rows, 1);
    out = nan(nGroups, width(T));

    for g = 1:nGroups
        idx = rows(g, :);
        idx = idx(isfinite(idx) & idx > 0);
        idx = idx(:);
        X = T{idx, :};
        out(g, :) = mean(X, 1, 'omitnan');
        idIdx = find(varNames == "POINTID", 1);
        if ~isempty(idIdx)
            out(g, idIdx) = T.POINTID(idx(1));
        end
    end

    pointTable = array2table(out, 'VariableNames', T.Properties.VariableNames);
end

function D = localComplete(T, vars)
    D = T(:, vars);
    D = D(all(isfinite(D{:,:}), 2), :);
end

function metrics = localCrossValidatedRF(T, predictorNames, targetName, useLog1pTarget)
    X = T{:, predictorNames};
    y = T.(targetName);
    [X, y] = localTransform(X, y, predictorNames);

    cv = cvpartition(size(X, 1), 'KFold', 5);
    yhat = nan(size(y));

    for k = 1:cv.NumTestSets
        tr = training(cv, k);
        te = test(cv, k);
        mu = mean(X(tr, :), 1);
        sig = std(X(tr, :), 0, 1);
        sig(sig == 0) = 1;
        Xtr = (X(tr, :) - mu) ./ sig;
        Xte = (X(te, :) - mu) ./ sig;

        if useLog1pTarget
            yTrain = log1p(y(tr));
        else
            yTrain = log10(max(y(tr), realmin));
        end

        rf = TreeBagger(300, Xtr, yTrain, ...
            'Method', 'regression', ...
            'MinLeafSize', 4, ...
            'NumPredictorsToSample', max(1, round(sqrt(numel(predictorNames)))));
        pred = predict(rf, Xte);
        if useLog1pTarget
            yhat(te) = exp(pred) - 1;
        else
            yhat(te) = 10 .^ pred;
        end
    end

    [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat);
    metrics = struct('RMSE', rmse, 'MAE', mae, 'R2_log10', r2, 'RSS_log10', rss, 'TSS_log10', tss);
end

function [X, y] = localTransform(X, y, predictorNames)
    skewed = ["OC","CaCO3","N","P","K","Coarse","Clay","CEC"];
    for j = 1:numel(predictorNames)
        if any(predictorNames(j) == skewed)
            X(:, j) = log1p(max(X(:, j), 0));
        end
    end
    y = max(y, realmin);
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
