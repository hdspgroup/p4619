function soil_ec_analysis()
% Analysis for non-linear soil EC modeling from data.mat
% Saves figures and a results MAT/text summary in the current folder.

    close all force;
    clc;

    paths = codex_paths('soil_ec_analysis');
    inFile = fullfile(paths.dataDir, 'data.mat');
    outDir = paths.analysisDir;
    s = load(inFile);
    data = s.data;
    varNames = matlab.lang.makeValidName(string(s.varNames), 'ReplacementStyle', 'delete');

    T = array2table(data, 'VariableNames', cellstr(varNames));

    predictorNames = ["pH_CaCl2","pH_H2O","OC","P","N","K","Coarse","Clay","Sand","Silt"];
    contextNames = ["sigma_w","T","H_2OVol","saturation"];
    targetNames = ["EC_sat","EC_dry"];

    skewedPredictors = ["OC","P","N","K","Coarse","Clay"];
    skewedTargets = ["EC_sat","EC_dry"];

    transformed = T;
    transformMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

    for name = skewedPredictors
        transformed.(name) = localLog1p(T.(name));
        transformMap(char(name)) = 'log1p';
    end

    for name = skewedTargets
        transformed.(name) = localLog1p(T.(name));
        transformMap(char(name)) = 'log1p';
    end

    % Summaries used later in the report.
    nRows = height(T);
    nMissing = varfun(@(x) sum(isnan(x)), T, 'OutputFormat', 'uniform');
    nUnique = zeros(1, width(T));
    for j = 1:width(T)
        x = T{:, j};
        nUnique(j) = numel(unique(x(~isnan(x))));
    end

    summaryTbl = table(string(T.Properties.VariableNames)', nMissing', nUnique', ...
        'VariableNames', {'Variable','Missing','Unique'});

    writetable(summaryTbl, fullfile(outDir, 'variable_summary.csv'));

    % Correlation heatmap on informative transformed variables.
    heatmapVars = ["EC_sat","EC_dry","pH_CaCl2","pH_H2O","OC","P","N","K","Coarse","Clay","Sand","Silt"];
    H = transformed(:, heatmapVars);
    corrMat = corr(table2array(H), 'Rows', 'pairwise', 'Type', 'Spearman');

    f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 800]);
    imagesc(corrMat);
    axis tight;
    axis equal;
    colormap(turbo(256));
    cb = colorbar;
    cb.Label.String = 'Spearman correlation';
    xticks(1:numel(heatmapVars));
    yticks(1:numel(heatmapVars));
    xticklabels(heatmapVars);
    yticklabels(heatmapVars);
    xtickangle(45);
    title('Correlation heatmap (transformed variables)');
    for r = 1:size(corrMat, 1)
        for c = 1:size(corrMat, 2)
            text(c, r, sprintf('%.2f', corrMat(r, c)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'k');
        end
    end
    exportgraphics(f1, fullfile(outDir, 'correlation_heatmap.png'), 'Resolution', 200);
    close(f1);

    % Scatter plot matrix. Moisture is constant in this dataset, so we omit it.
    scatterVars = ["EC_sat","EC_dry","pH_H2O","pH_CaCl2","OC","N","P","K","Clay","Sand","Silt"];
    S = transformed(:, scatterVars);
    keepScatter = all(~ismissing(S), 2);
    S = S(keepScatter, :);
    nScatter = min(height(S), 2500);
    rng(7);
    if height(S) > nScatter
        S = S(randperm(height(S), nScatter), :);
    end
    f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1200]);
    plotmatrix(table2array(S));
    sgtitle('Scatter plot matrix (transformed variables)');
    exportgraphics(f2, fullfile(outDir, 'scatter_matrix.png'), 'Resolution', 180);
    close(f2);

    % Model comparison for both targets.
    results = struct();
    for targetName = targetNames
        [metrics, predTbl, importanceTbl, sensitivityTbl, physicsCoefTbl] = localCompareModels(transformed, predictorNames, targetName);
        results.(targetName).metrics = metrics;
        results.(targetName).predictions = predTbl;
        results.(targetName).importance = importanceTbl;
        results.(targetName).sensitivity = sensitivityTbl;
        results.(targetName).physicsCoefficients = physicsCoefTbl;

        writetable(predTbl, fullfile(outDir, sprintf('predictions_%s.csv', targetName)));
        writetable(importanceTbl, fullfile(outDir, sprintf('rf_importance_%s.csv', targetName)));
        writetable(sensitivityTbl, fullfile(outDir, sprintf('sensitivity_%s.csv', targetName)));
        writetable(physicsCoefTbl, fullfile(outDir, sprintf('physics_coefficients_%s.csv', targetName)));

        fSens = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 500]);
        bar(categorical(sensitivityTbl.Predictor), sensitivityTbl.TSSScaledUncertaintyContribution);
        ylabel('Contribution / TSS');
        title(sprintf('Gradient-based uncertainty contribution for %s', targetName));
        grid on;
        exportgraphics(fSens, fullfile(outDir, sprintf('sensitivity_%s.png', targetName)), 'Resolution', 180);
        close(fSens);

        localSavePredictedVsActualPlot(predTbl, metrics, targetName, outDir);
    end

    save(fullfile(outDir, 'soil_ec_results.mat'), 'results', 'summaryTbl', 'corrMat', 'heatmapVars');

    fid = fopen(fullfile(outDir, 'model_summary.txt'), 'w');
    fprintf(fid, 'Soil EC modeling summary\n');
    fprintf(fid, 'Rows in file: %d\n\n', nRows);
    fprintf(fid, 'Constant-context variables:\n');
    for name = contextNames
        fprintf(fid, '  %s unique count = %d\n', name, summaryTbl.Unique(summaryTbl.Variable == name));
    end
    fprintf(fid, '\nTransformation rules:\n');
    keys = string(transformMap.keys);
    for i = 1:numel(keys)
        fprintf(fid, '  %s -> %s\n', keys(i), transformMap(char(keys(i))));
    end
    fprintf(fid, '\nModel comparison (model selection and variance metrics in log10 space; RMSE/MAE in original scale):\n');

    for targetName = targetNames
        m = results.(targetName).metrics;
        fprintf(fid, '\nTarget: %s\n', targetName);
        fprintf(fid, '  Complete cases used: %d\n', m.N);
        fprintf(fid, '  Total sum of squares (TSS, log10 space) = %.6g\n', m.TSS);
        fprintf(fid, '  Physics-informed: RMSE = %.6g, MAE = %.6g, RSS = %.6g, ESS = %.6g, R^2 = %.4f\n', ...
            m.Physics_RMSE, m.Physics_MAE, m.Physics_RSS, m.Physics_ESS, m.Physics_R2);
        fprintf(fid, '  Power law: RMSE = %.6g, MAE = %.6g, RSS = %.6g, ESS = %.6g, R^2 = %.4f\n', ...
            m.Power_RMSE, m.Power_MAE, m.Power_RSS, m.Power_ESS, m.Power_R2);
        fprintf(fid, '  Poly degree 2: RMSE = %.6g, MAE = %.6g, RSS = %.6g, ESS = %.6g, R^2 = %.4f\n', ...
            m.Poly2_RMSE, m.Poly2_MAE, m.Poly2_RSS, m.Poly2_ESS, m.Poly2_R2);
        fprintf(fid, '  Poly degree 3: RMSE = %.6g, MAE = %.6g, RSS = %.6g, ESS = %.6g, R^2 = %.4f\n', ...
            m.Poly3_RMSE, m.Poly3_MAE, m.Poly3_RSS, m.Poly3_ESS, m.Poly3_R2);
        fprintf(fid, '  Random forest: RMSE = %.6g, MAE = %.6g, RSS = %.6g, ESS = %.6g, R^2 = %.4f\n', ...
            m.RF_RMSE, m.RF_MAE, m.RF_RSS, m.RF_ESS, m.RF_R2);
        fprintf(fid, '  Best polynomial degree: %d\n', m.BestPolynomialDegree);
        fprintf(fid, '  Best overall model: %s\n', m.BestModel);
        fprintf(fid, '  Highest TSS-scaled uncertainty contributors in best model:\n');
        sens = results.(targetName).sensitivity;
        topN = min(5, height(sens));
        for i = 1:topN
            fprintf(fid, '    %s: contribution/TSS = %.6g, mean|dY/dX| = %.6g, uncertainty proxy = %.6g\n', ...
                sens.Predictor(i), sens.TSSScaledUncertaintyContribution(i), ...
                sens.MeanAbsGradient(i), sens.UncertaintyProxy(i));
        end
    end
    fclose(fid);

    disp('Analysis complete. Outputs written to current folder.');
end

function [metrics, predTbl, importanceTbl, sensitivityTbl, physicsCoefTbl] = localCompareModels(T, predictorNames, targetName)
    useVars = [predictorNames, targetName];
    D = T(:, useVars);
    keep = all(~ismissing(D), 2);
    D = D(keep, :);

    X = D{:, predictorNames};
    y = D{:, targetName};
    yOriginal = localInvertTargetTransform(y, targetName);
    yLog10 = log10(max(yOriginal, realmin));

    cv = cvpartition(size(X, 1), 'KFold', 5);
    yhatPhysics = nan(size(y));
    yhatPower = nan(size(y));
    yhat2 = nan(size(y));
    yhat3 = nan(size(y));
    yhatRF = nan(size(y));
    rfImportance = zeros(cv.NumTestSets, numel(predictorNames));
    physicsCoefMat = [];
    physicsFeatureNames = string.empty;

    for k = 1:cv.NumTestSets
        trainIdx = training(cv, k);
        testIdx = test(cv, k);

        mu = mean(X(trainIdx, :), 1);
        sigma = std(X(trainIdx, :), 0, 1);
        sigma(sigma == 0) = 1;

        XTrain = (X(trainIdx, :) - mu) ./ sigma;
        XTest = (X(testIdx, :) - mu) ./ sigma;
        yTrain = y(trainIdx);
        yTrainOriginal = yOriginal(trainIdx);
        yTrainLog10 = log10(max(yTrainOriginal, realmin));

        powerTblTrain = array2table(XTrain);
        powerTblTest = array2table(XTest);
        for j = 1:width(powerTblTrain)
            vName = sprintf('x%d', j);
            powerTblTrain.Properties.VariableNames{j} = vName;
            powerTblTest.Properties.VariableNames{j} = vName;
        end

        mdlPower = fitlm(powerTblTrain, yTrainLog10);
        [XTrainPhys, physicsFeatureNames] = localPhysicsFeatures(X(trainIdx, :), predictorNames);
        XTestPhys = localPhysicsFeatures(X(testIdx, :), predictorNames);
        [BPhys, fitInfoPhys] = lasso(XTrainPhys, yTrainLog10, ...
            'CV', 5, 'Alpha', 0.6, 'Standardize', true);
        idxPhys = fitInfoPhys.IndexMinMSE;
        coefPhys = BPhys(:, idxPhys);
        interceptPhys = fitInfoPhys.Intercept(idxPhys);

        mdl2 = fitlm(localPolyFeatures(XTrain, 2), yTrain);
        mdl3 = fitlm(localPolyFeatures(XTrain, 3), yTrain);

        predPhysicsLog10 = XTestPhys * coefPhys + interceptPhys;
        predPowerLog10 = predict(mdlPower, powerTblTest);
        pred2 = predict(mdl2, localPolyFeatures(XTest, 2));
        pred3 = predict(mdl3, localPolyFeatures(XTest, 3));

        rf = TreeBagger(300, XTrain, yTrain, ...
            'Method', 'regression', ...
            'MinLeafSize', 5, ...
            'NumPredictorsToSample', max(1, round(sqrt(numel(predictorNames)))), ...
            'OOBPredictorImportance', 'on');
        predRF = predict(rf, XTest);

        yhatPhysics(testIdx) = 10.^predPhysicsLog10;
        yhatPower(testIdx) = 10.^predPowerLog10;
        yhat2(testIdx) = pred2;
        yhat3(testIdx) = pred3;
        yhatRF(testIdx) = predRF;
        rfImportance(k, :) = rf.OOBPermutedPredictorDeltaError;
        physicsCoefMat = [physicsCoefMat, [interceptPhys; coefPhys]]; %#ok<AGROW>
    end

    yhatPhysicsOriginal = yhatPhysics;
    yhatPowerOriginal = yhatPower;
    yhat2Original = localInvertTargetTransform(yhat2, targetName);
    yhat3Original = localInvertTargetTransform(yhat3, targetName);
    yhatRFOriginal = localInvertTargetTransform(yhatRF, targetName);

    metrics = struct();
    metrics.N = numel(yOriginal);
    [metrics.Physics_RMSE, metrics.Physics_MAE, metrics.Physics_R2, metrics.Physics_RSS, metrics.TSS] = localRegressionMetrics(yOriginal, yhatPhysicsOriginal);
    metrics.Physics_ESS = metrics.TSS - metrics.Physics_RSS;
    [metrics.Power_RMSE, metrics.Power_MAE, metrics.Power_R2, metrics.Power_RSS, metrics.TSS] = localRegressionMetrics(yOriginal, yhatPowerOriginal);
    metrics.Power_ESS = metrics.TSS - metrics.Power_RSS;
    [metrics.Poly2_RMSE, metrics.Poly2_MAE, metrics.Poly2_R2, metrics.Poly2_RSS] = localRegressionMetrics(yOriginal, yhat2Original);
    metrics.Poly2_ESS = metrics.TSS - metrics.Poly2_RSS;
    [metrics.Poly3_RMSE, metrics.Poly3_MAE, metrics.Poly3_R2, metrics.Poly3_RSS] = localRegressionMetrics(yOriginal, yhat3Original);
    metrics.Poly3_ESS = metrics.TSS - metrics.Poly3_RSS;
    [metrics.RF_RMSE, metrics.RF_MAE, metrics.RF_R2, metrics.RF_RSS] = localRegressionMetrics(yOriginal, yhatRFOriginal);
    metrics.RF_ESS = metrics.TSS - metrics.RF_RSS;

    if metrics.Poly3_R2 >= metrics.Poly2_R2
        metrics.BestPolynomialDegree = 3;
        bestPolyR2 = metrics.Poly3_R2;
    else
        metrics.BestPolynomialDegree = 2;
        bestPolyR2 = metrics.Poly2_R2;
    end

    [~, bestIdx] = max([metrics.Physics_R2, metrics.Power_R2, bestPolyR2, metrics.RF_R2]);
    if bestIdx == 1
        metrics.BestModel = 'Physics-informed';
    elseif bestIdx == 2
        metrics.BestModel = 'Power law';
    elseif bestIdx == 3
        metrics.BestModel = sprintf('Polynomial degree %d', metrics.BestPolynomialDegree);
    else
        metrics.BestModel = 'Random Forest';
    end

    predTbl = table(yOriginal, yhatPhysicsOriginal, yhatPowerOriginal, yhat2Original, yhat3Original, yhatRFOriginal, ...
        'VariableNames', {'Observed','PhysicsInformed','PowerLaw','Poly2','Poly3','RandomForest'});

    importanceTbl = table(predictorNames', mean(rfImportance, 1)', std(rfImportance, 0, 1)', ...
        'VariableNames', {'Predictor','MeanOOBImportance','StdOOBImportance'});
    importanceTbl = sortrows(importanceTbl, 'MeanOOBImportance', 'descend');

    sensitivityTbl = localSensitivityAnalysis(X, yOriginal, predictorNames, targetName, metrics);
    physicsCoefTbl = table(["Intercept"; physicsFeatureNames(:)], mean(physicsCoefMat, 2), std(physicsCoefMat, 0, 2), ...
        'VariableNames', {'Term','MeanCoefficient','StdCoefficient'});
    physicsCoefTbl = sortrows(physicsCoefTbl(2:end, :), 'MeanCoefficient', 'descend');
    physicsCoefTbl = [table("Intercept", mean(physicsCoefMat(1, :)), std(physicsCoefMat(1, :), 0, 2), ...
        'VariableNames', {'Term','MeanCoefficient','StdCoefficient'}); physicsCoefTbl];
end

function Z = localPolyFeatures(X, degree)
    n = size(X, 1);
    p = size(X, 2);
    Z = [X, X.^2];

    for i = 1:p - 1
        for j = i + 1:p
            Z = [Z, X(:, i) .* X(:, j)]; %#ok<AGROW>
        end
    end

    if degree >= 3
        Z = [Z, X.^3];
        for i = 1:p
            for j = 1:p
                if i ~= j
                    Z = [Z, (X(:, i).^2) .* X(:, j)]; %#ok<AGROW>
                end
            end
        end
    end

    Z = array2table(Z);
    for j = 1:width(Z)
        Z.Properties.VariableNames{j} = sprintf('x%d', j);
    end
end

function x = localLog1p(x)
    x = log1p(x);
end

function x = localInvertTargetTransform(x, targetName)
    if any(strcmp(targetName, ["EC_sat","EC_dry"]))
        x = exp(x) - 1;
        x(~isfinite(x)) = NaN;
    end
end

function [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat)
    keep = isfinite(y) & isfinite(yhat);
    y = y(keep);
    yhat = yhat(keep);
    rmse = sqrt(mean((y - yhat).^2));
    mae = mean(abs(y - yhat));
    ySafe = max(y, realmin);
    yhatSafe = max(yhat, realmin);
    yLog = log10(ySafe);
    yhatLog = log10(yhatSafe);
    tss = sum((yLog - mean(yLog)).^2);
    rss = sum((yLog - yhatLog).^2);
    r2 = 1 - rss / tss;
end

function sensitivityTbl = localSensitivityAnalysis(X, yOriginal, predictorNames, targetName, metrics)
    mu = mean(X, 1);
    sigma = std(X, 0, 1);
    sigma(sigma == 0) = 1;
    XScaled = (X - mu) ./ sigma;

    if strcmp(metrics.BestModel, 'Random Forest')
        finalModel = TreeBagger(300, XScaled, log1p(yOriginal), ...
            'Method', 'regression', ...
            'MinLeafSize', 5, ...
            'NumPredictorsToSample', max(1, round(sqrt(numel(predictorNames)))));
        predictFcn = @(Xin) exp(predict(finalModel, (Xin - mu) ./ sigma)) - 1;
    elseif strcmp(metrics.BestModel, 'Physics-informed')
        [XPhys, physNames] = localPhysicsFeatures(X, predictorNames);
        [BPhys, fitInfoPhys] = lasso(XPhys, log10(max(yOriginal, realmin)), ...
            'CV', 5, 'Alpha', 0.6, 'Standardize', true);
        idxPhys = fitInfoPhys.IndexMinMSE;
        coefPhys = BPhys(:, idxPhys);
        interceptPhys = fitInfoPhys.Intercept(idxPhys);
        predictFcn = @(Xin) 10.^(localPhysicsFeatures(Xin, predictorNames) * coefPhys + interceptPhys);
    elseif strcmp(metrics.BestModel, 'Power law')
        powerTbl = array2table(XScaled);
        for j = 1:width(powerTbl)
            powerTbl.Properties.VariableNames{j} = sprintf('x%d', j);
        end
        finalModel = fitlm(powerTbl, log10(max(yOriginal, realmin)));
        predictFcn = @(Xin) 10.^predict(finalModel, localNamedTable((Xin - mu) ./ sigma));
    else
        degree = metrics.BestPolynomialDegree;
        finalModel = fitlm(localPolyFeatures(XScaled, degree), log1p(yOriginal));
        predictFcn = @(Xin) exp(predict(finalModel, localPolyFeatures((Xin - mu) ./ sigma, degree))) - 1;
    end

    p = size(X, 2);
    meanAbsGrad = zeros(p, 1);
    medianAbsGrad = zeros(p, 1);
    p95AbsGrad = zeros(p, 1);
    uncertaintyProxy = zeros(p, 1);
    tssContribution = zeros(p, 1);

    for j = 1:p
        xj = X(:, j);
        step = max(1e-6, 0.01 * iqr(xj));
        if step == 0
            step = max(1e-6, 0.01 * std(xj));
        end
        Xplus = X;
        Xminus = X;
        Xplus(:, j) = Xplus(:, j) + step;
        Xminus(:, j) = max(0, Xminus(:, j) - step);

        yPlus = predictFcn(Xplus);
        yMinus = predictFcn(Xminus);
        grad = (yPlus - yMinus) ./ (Xplus(:, j) - Xminus(:, j));
        absGrad = abs(grad);

        meanAbsGrad(j) = mean(absGrad, 'omitnan');
        medianAbsGrad(j) = median(absGrad, 'omitnan');
        p95AbsGrad(j) = prctile(absGrad, 95);
        uncertaintyProxy(j) = 1.4826 * mad(xj, 1);
        tssContribution(j) = mean((grad .* uncertaintyProxy(j)).^2, 'omitnan');
    end

    sensitivityTbl = table(predictorNames', meanAbsGrad, medianAbsGrad, p95AbsGrad, ...
        uncertaintyProxy, tssContribution, tssContribution ./ metrics.TSS, ...
        'VariableNames', {'Predictor','MeanAbsGradient','MedianAbsGradient', ...
        'P95AbsGradient','UncertaintyProxy','UncertaintyContribution','TSSScaledUncertaintyContribution'});
    sensitivityTbl = sortrows(sensitivityTbl, 'TSSScaledUncertaintyContribution', 'descend');
end

function localSavePredictedVsActualPlot(predTbl, metrics, targetName, outDir)
    if strcmp(metrics.BestModel, 'Random Forest')
        bestPred = predTbl.RandomForest;
        modelLabel = 'Random Forest';
    elseif strcmp(metrics.BestModel, 'Physics-informed')
        bestPred = predTbl.PhysicsInformed;
        modelLabel = 'Physics-informed';
    elseif strcmp(metrics.BestModel, 'Power law')
        bestPred = predTbl.PowerLaw;
        modelLabel = 'Power law';
    else
        degree = metrics.BestPolynomialDegree;
        if degree == 2
            bestPred = predTbl.Poly2;
        else
            bestPred = predTbl.Poly3;
        end
        modelLabel = sprintf('Polynomial degree %d', degree);
    end

    observed = predTbl.Observed;
    keep = isfinite(observed) & isfinite(bestPred);
    observed = observed(keep);
    bestPred = bestPred(keep);

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 700]);
    minPositive = min([observed(observed > 0); bestPred(bestPred > 0)]);
    observed(observed <= 0) = minPositive * 0.5;
    bestPred(bestPred <= 0) = minPositive * 0.5;

    scatter(observed, bestPred, 12, 'filled', 'MarkerFaceAlpha', 0.28);
    hold on;
    minVal = min([observed; bestPred]);
    maxVal = max([observed; bestPred]);
    plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Actual');
    ylabel('Predicted');
    title(sprintf('Predicted vs Actual for %s (%s, log-log axes)', targetName, modelLabel));
    txt = sprintf('R^2 = %.3f | RMSE = %.3g | MAE = %.3g', ...
        localGetMetric(metrics, modelLabel, 'R2'), ...
        localGetMetric(metrics, modelLabel, 'RMSE'), ...
        localGetMetric(metrics, modelLabel, 'MAE'));
    text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 6);
    grid on;
    box on;
    exportgraphics(f, fullfile(outDir, sprintf('predicted_vs_actual_%s.png', targetName)), 'Resolution', 200);
    close(f);
end

function value = localGetMetric(metrics, modelLabel, metricName)
    if strcmp(modelLabel, 'Random Forest')
        prefix = 'RF';
    elseif strcmp(modelLabel, 'Physics-informed')
        prefix = 'Physics';
    elseif strcmp(modelLabel, 'Power law')
        prefix = 'Power';
    elseif contains(modelLabel, 'degree 2')
        prefix = 'Poly2';
    else
        prefix = 'Poly3';
    end
    value = metrics.([prefix '_' metricName]);
end

function T = localNamedTable(X)
    T = array2table(X);
    for j = 1:width(T)
        T.Properties.VariableNames{j} = sprintf('x%d', j);
    end
end

function [Z, featureNames] = localPhysicsFeatures(X, predictorNames)
    idx = @(name) find(strcmp(predictorNames, name), 1);

    pHCa = X(:, idx("pH_CaCl2"));
    pHH2O = X(:, idx("pH_H2O"));
    OC = X(:, idx("OC"));
    P = X(:, idx("P"));
    N = X(:, idx("N"));
    K = X(:, idx("K"));
    Coarse = X(:, idx("Coarse"));
    Clay = X(:, idx("Clay"));
    Sand = X(:, idx("Sand"));
    Silt = X(:, idx("Silt"));

    lOC = log1p(OC);
    lP = log1p(P);
    lN = log1p(N);
    lK = log1p(K);
    lCoarse = log1p(Coarse);
    lClay = log1p(Clay);
    lSand = log1p(Sand);
    lSilt = log1p(Silt);
    lFine = log1p(Clay + Silt);
    lTextureRatio = log((Clay + Silt + 1) ./ (Sand + Coarse + 1));
    deltaPH = pHH2O - pHCa;
    meanPH = 0.5 * (pHH2O + pHCa);
    finesShare = (Clay + Silt) ./ max(Clay + Silt + Sand + Coarse, 1);
    coarseShare = (Sand + Coarse) ./ max(Clay + Silt + Sand + Coarse, 1);

    % Physics-informed proxies:
    % 1) dissolved electrolyte content in extract
    % 2) surface conduction driven by fine particles and organic matter
    % 3) dilution/blocking by coarse texture
    Z = [ ...
        pHCa, pHH2O, meanPH, deltaPH, ...
        lOC, lP, lN, lK, lClay, lSilt, lSand, lCoarse, lFine, lTextureRatio, finesShare, coarseShare, ...
        lClay .* lK, ...
        lFine .* lK, ...
        lClay .* lP, ...
        lClay .* lOC, lFine .* lOC, ...
        lSand .* lK, lCoarse .* lK, ...
        meanPH .* lK, meanPH .* lP, ...
        deltaPH .* lClay, deltaPH .* lFine, ...
        pHH2O .* lOC, pHCa .* lClay, ...
        lClay .* lN];

    featureNames = [ ...
        "pH_CaCl2","pH_H2O","mean_pH","delta_pH", ...
        "log1p_OC","log1p_P","log1p_N","log1p_K","log1p_Clay","log1p_Silt","log1p_Sand","log1p_Coarse", ...
        "log1p_Fine","log_TextureRatio","FineShare","CoarseShare", ...
        "Clay_x_K","Fine_x_K","Clay_x_P", ...
        "Clay_x_OC","Fine_x_OC", ...
        "Sand_x_K","Coarse_x_K", ...
        "meanpH_x_K","meanpH_x_P", ...
        "deltapH_x_Clay","deltapH_x_Fine","pHH2O_x_OC","pHCa_x_Clay", ...
        "Clay_x_TotalN"];
end
