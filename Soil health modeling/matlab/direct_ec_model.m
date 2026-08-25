function direct_ec_model()
% Direct predictive modeling for saturated soil EC using current variables.
% Metrics for model selection are computed in log10 space.

    close all force;
    clc;

    paths = codex_paths('direct_ec');
    outDir = paths.analysisDir;

    s = load(fullfile(paths.dataDir, 'data.mat'));
    rawNames = string(s.varNames);
    names = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
    T = array2table(s.data, 'VariableNames', cellstr(names));

    predictorNames = ["pH_CaCl2","pH_H2O","OC","P","N","K","Coarse","Clay","Sand","Silt"];
    targetName = "EC_sat";

    skewedPredictors = ["OC","P","N","K","Coarse","Clay"];
    for v = skewedPredictors
        T.(v) = log1p(T.(v));
    end

    D = T(:, [predictorNames, targetName]);
    keep = all(~ismissing(D), 2);
    D = D(keep, :);

    X = D{:, predictorNames};
    y = D{:, targetName};

    cv = cvpartition(size(X, 1), 'KFold', 5);
    yhatPoly2 = nan(size(y));
    yhatRF = nan(size(y));
    yhatBoost = nan(size(y));

    for k = 1:cv.NumTestSets
        trainIdx = training(cv, k);
        testIdx = test(cv, k);

        XTrain = X(trainIdx, :);
        XTest = X(testIdx, :);
        yTrain = y(trainIdx);

        mu = mean(XTrain, 1);
        sigma = std(XTrain, 0, 1);
        sigma(sigma == 0) = 1;
        XTrainScaled = (XTrain - mu) ./ sigma;
        XTestScaled = (XTest - mu) ./ sigma;

        mdlPoly2 = fitlm(localPolyFeatures(XTrainScaled, 2), log1p(yTrain));
        predPoly2 = exp(predict(mdlPoly2, localPolyFeatures(XTestScaled, 2))) - 1;

        mdlRF = TreeBagger(400, XTrainScaled, log1p(yTrain), ...
            'Method', 'regression', ...
            'MinLeafSize', 4, ...
            'NumPredictorsToSample', max(1, round(sqrt(numel(predictorNames)))));
        predRF = exp(predict(mdlRF, XTestScaled)) - 1;

        tblTrain = array2table(XTrain, 'VariableNames', cellstr(predictorNames));
        tblTest = array2table(XTest, 'VariableNames', cellstr(predictorNames));
        mdlBoost = fitrensemble(tblTrain, log10(max(yTrain, realmin)), ...
            'Method', 'LSBoost', ...
            'NumLearningCycles', 300, ...
            'Learners', templateTree('MinLeafSize', 8, 'NumVariablesToSample', numel(predictorNames)));
        predBoost = 10.^predict(mdlBoost, tblTest);

        yhatPoly2(testIdx) = predPoly2;
        yhatRF(testIdx) = predRF;
        yhatBoost(testIdx) = predBoost;
    end

    metrics = struct();
    [metrics.Poly2_RMSE, metrics.Poly2_MAE, metrics.Poly2_R2, metrics.Poly2_RSS, metrics.TSS] = localRegressionMetrics(y, yhatPoly2);
    [metrics.RF_RMSE, metrics.RF_MAE, metrics.RF_R2, metrics.RF_RSS] = localRegressionMetrics(y, yhatRF);
    [metrics.Boost_RMSE, metrics.Boost_MAE, metrics.Boost_R2, metrics.Boost_RSS] = localRegressionMetrics(y, yhatBoost);

    [~, bestIdx] = max([metrics.Poly2_R2, metrics.RF_R2, metrics.Boost_R2]);
    labels = ["Polynomial degree 2","Random Forest","LSBoost"];
    metrics.BestModel = labels(bestIdx);

    predTbl = table(y, yhatPoly2, yhatRF, yhatBoost, ...
        'VariableNames', {'Observed','Poly2','RandomForest','LSBoost'});
    writetable(predTbl, fullfile(outDir, 'direct_ec_predictions.csv'));

    summaryTbl = table( ...
        ["Polynomial degree 2"; "Random Forest"; "LSBoost"], ...
        [metrics.Poly2_RMSE; metrics.RF_RMSE; metrics.Boost_RMSE], ...
        [metrics.Poly2_MAE; metrics.RF_MAE; metrics.Boost_MAE], ...
        [metrics.Poly2_RSS; metrics.RF_RSS; metrics.Boost_RSS], ...
        [metrics.Poly2_R2; metrics.RF_R2; metrics.Boost_R2], ...
        'VariableNames', {'Model','RMSE','MAE','RSS_log10','R2_log10'});
    writetable(summaryTbl, fullfile(outDir, 'direct_ec_summary.csv'));

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 820 700]);
    switch metrics.BestModel
        case "Polynomial degree 2"
            pred = yhatPoly2;
        case "Random Forest"
            pred = yhatRF;
        otherwise
            pred = yhatBoost;
    end
    obs = max(y, realmin);
    pred = max(pred, realmin);
    scatter(obs, pred, 12, 'filled', 'MarkerFaceAlpha', 0.28);
    hold on;
    lo = min([obs; pred]);
    hi = max([obs; pred]);
    plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Actual EC_{sat}');
    ylabel('Predicted EC_{sat}');
    title(sprintf('Direct EC model: %s', metrics.BestModel));
    txt = sprintf('R^2_{log10} = %.3f | RMSE = %.3g | MAE = %.3g', ...
        max([metrics.Poly2_R2, metrics.RF_R2, metrics.Boost_R2]), ...
        summaryTbl.RMSE(bestIdx), summaryTbl.MAE(bestIdx));
    text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 6);
    grid on;
    exportgraphics(f, fullfile(outDir, 'direct_ec_predicted_vs_actual.png'), 'Resolution', 200);
    close(f);

    fid = fopen(fullfile(outDir, 'direct_ec_summary.txt'), 'w');
    fprintf(fid, 'Direct EC model summary\n');
    fprintf(fid, 'Complete cases: %d\n', numel(y));
    fprintf(fid, 'TSS in log10 space: %.6g\n\n', metrics.TSS);
    fprintf(fid, 'Polynomial degree 2: RMSE = %.6g, MAE = %.6g, RSS = %.6g, R^2 = %.4f\n', ...
        metrics.Poly2_RMSE, metrics.Poly2_MAE, metrics.Poly2_RSS, metrics.Poly2_R2);
    fprintf(fid, 'Random Forest: RMSE = %.6g, MAE = %.6g, RSS = %.6g, R^2 = %.4f\n', ...
        metrics.RF_RMSE, metrics.RF_MAE, metrics.RF_RSS, metrics.RF_R2);
    fprintf(fid, 'LSBoost: RMSE = %.6g, MAE = %.6g, RSS = %.6g, R^2 = %.4f\n', ...
        metrics.Boost_RMSE, metrics.Boost_MAE, metrics.Boost_RSS, metrics.Boost_R2);
    fprintf(fid, '\nBest direct model: %s\n', metrics.BestModel);
    fclose(fid);

    % Fit final degree-2 polynomial on all complete cases and export the explicit formula.
    mu = mean(X, 1);
    sigma = std(X, 0, 1);
    sigma(sigma == 0) = 1;
    XScaled = (X - mu) ./ sigma;
    polyTerms = localPolyFeatures(XScaled, 2);
    finalPoly = fitlm(polyTerms, log1p(y));

    predictorTbl = table(predictorNames', mu', sigma', ...
        'VariableNames', {'Predictor','Mean','StdDev'});
    writetable(predictorTbl, fullfile(outDir, 'direct_ec_poly2_scaling.csv'));

    coeffTbl = finalPoly.Coefficients;
    coeffNames = ["Intercept"; string(polyTerms.Properties.VariableNames(:))];
    coeffOut = table(coeffNames, coeffTbl.Estimate, coeffTbl.SE, coeffTbl.tStat, coeffTbl.pValue, ...
        'VariableNames', {'Term','Estimate','SE','tStat','pValue'});
    writetable(coeffOut, fullfile(outDir, 'direct_ec_poly2_coefficients.csv'));

    formulaText = localBuildFormulaText(predictorNames, mu, sigma, coeffOut);
    fid = fopen(fullfile(outDir, 'direct_ec_poly2_formula.txt'), 'w');
    fprintf(fid, '%s', formulaText);
    fclose(fid);

    save(fullfile(outDir, 'direct_ec_model_results.mat'), 'metrics', 'summaryTbl');
    disp('Direct EC modeling complete.');
end

function Z = localPolyFeatures(X, degree)
    p = size(X, 2);
    Z = [X, X.^2];
    for i = 1:p - 1
        for j = i + 1:p
            Z = [Z, X(:, i) .* X(:, j)]; %#ok<AGROW>
        end
    end
    if degree >= 3
        Z = [Z, X.^3]; %#ok<UNRCH>
    end
    Z = array2table(Z);
    for j = 1:width(Z)
        Z.Properties.VariableNames{j} = sprintf('x%d', j);
    end
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

function txt = localBuildFormulaText(predictorNames, mu, sigma, coeffOut)
    lines = strings(0,1);
    lines(end+1) = "Direct EC degree-2 polynomial formula";
    lines(end+1) = "";
    lines(end+1) = "1. Transform raw predictors:";
    lines(end+1) = "   z1 = pH_CaCl2";
    lines(end+1) = "   z2 = pH_H2O";
    lines(end+1) = "   z3 = log(1 + OC)";
    lines(end+1) = "   z4 = log(1 + P)";
    lines(end+1) = "   z5 = log(1 + N)";
    lines(end+1) = "   z6 = log(1 + K)";
    lines(end+1) = "   z7 = log(1 + Coarse)";
    lines(end+1) = "   z8 = log(1 + Clay)";
    lines(end+1) = "   z9 = Sand";
    lines(end+1) = "   z10 = Silt";
    lines(end+1) = "";
    lines(end+1) = "2. Standardize predictors:";
    for i = 1:numel(predictorNames)
        lines(end+1) = sprintf("   x%d = (z%d - %.15g) / %.15g", i, i, mu(i), sigma(i));
    end
    lines(end+1) = "";
    lines(end+1) = "3. Compute u_hat = log(1 + EC_sat_hat):";

    expr = "u_hat = ";
    for i = 1:height(coeffOut)
        term = coeffOut.Term(i);
        beta = coeffOut.Estimate(i);
        if term == "Intercept"
            expr = expr + sprintf("%.15g", beta);
        else
            matlabTerm = localReadableTerm(term, predictorNames);
            if beta >= 0
                expr = expr + sprintf(" + %.15g*(%s)", beta, matlabTerm);
            else
                expr = expr + sprintf(" - %.15g*(%s)", abs(beta), matlabTerm);
            end
        end
    end
    lines(end+1) = "   " + expr;
    lines(end+1) = "";
    lines(end+1) = "4. Back-transform:";
    lines(end+1) = "   EC_sat_hat = exp(u_hat) - 1";
    txt = strjoin(lines, newline);
end

function out = localReadableTerm(term, predictorNames)
    if ismissing(term) || strlength(term) == 0
        out = "unknown";
        return;
    end
    if term == "Intercept"
        out = "1";
        return;
    end
    % Terms follow the column order created by localPolyFeatures:
    names = "x" + string(1:numel(predictorNames));
    p = numel(predictorNames);
    terms = strings(0,1);
    for i = 1:p
        terms(end+1) = names(i); %#ok<AGROW>
    end
    for i = 1:p
        terms(end+1) = names(i) + "^2"; %#ok<AGROW>
    end
    for i = 1:p-1
        for j = i+1:p
            terms(end+1) = names(i) + "*" + names(j); %#ok<AGROW>
        end
    end
    raw = char(extractAfter(term, "x"));
    pos = sscanf(raw, "%d");
    if ~isempty(pos) && pos >= 1 && pos <= numel(terms)
        out = char(terms(pos));
    else
        out = char(term);
    end
end
