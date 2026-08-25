function model = fit_inverse_soil_health_model(T, options)
% Fit an inverse soil-health model:
%   sensor variables + context variables -> soil-health state variables
%
% Sensor variables (measured in practice):
%   EC, pH, Moisture
%
% Context variables (soil modifiers / conditioning variables):
%   Clay, Silt, Sand, Coarse, CaCO3
%
% Soil-health state targets:
%   N, P, K, CEC, OC, Moisture
%
% The workflow supports missing optional variables. Moisture can be:
%   1) an input if it is measured
%   2) an output if it is not directly measured but available as a label
%   3) absent entirely
%
% A separate model is fit for each available target.
%
% Example:
%   T = readtable('my_inverse_training_data.csv');
%   model = fit_inverse_soil_health_model(T);

    arguments
        T table
        options.AnalysisName string = "inverse_soil_health"
        options.ECName string = "EC"
        options.ConvertECdSmToSm logical = false
    end

    paths = codex_paths(options.AnalysisName);
    outDir = paths.analysisDir;

    T = localNormalizeVariableNames(T);

    if ismember(options.ECName, string(T.Properties.VariableNames))
        ecName = options.ECName;
    elseif ismember("EC_sat", string(T.Properties.VariableNames))
        ecName = "EC_sat";
    elseif ismember("EC", string(T.Properties.VariableNames))
        ecName = "EC";
    else
        error('No supported EC input variable found.');
    end

    if options.ConvertECdSmToSm && ismember(ecName, string(T.Properties.VariableNames))
        T.(ecName) = T.(ecName) .* 0.1;
    end

    expectedTargets = ["N","P","K","CEC","OC","Moisture"];
    contextCandidates = ["Clay","Silt","Sand","Coarse","CaCO3"];
    sensorCandidates = [ecName, "pH_CaCl2", "pH_H2O", "Moisture"];

    contextVars = intersect(contextCandidates, string(T.Properties.VariableNames), 'stable');
    sensorVars = intersect(sensorCandidates, string(T.Properties.VariableNames), 'stable');
    availableTargets = intersect(expectedTargets, string(T.Properties.VariableNames), 'stable');

    % If moisture is used as a sensor input, do not also predict it.
    targetVars = setdiff(availableTargets, intersect(availableTargets, sensorVars), 'stable');

    if isempty(targetVars)
        error('No inverse-model targets are available after accounting for measured inputs.');
    end

    summaryRows = strings(0, 1);
    summaryN = [];
    summaryRMSE = [];
    summaryMAE = [];
    summaryR2 = [];

    model = struct();
    model.type = "inverse_soil_health";
    model.sensorVars = sensorVars;
    model.contextVars = contextVars;
    model.targetVars = targetVars;
    model.targetModels = struct();

    for t = targetVars
        predictors = [sensorVars, contextVars];
        D = T(:, [predictors, t]);
        D = D(all(~ismissing(D), 2), :);
        D = D(all(isfinite(D{:,:}), 2), :);

        if height(D) < 50
            continue;
        end

        [X, featureNames] = localInverseFeatures(D, sensorVars, contextVars);
        y = D.(t);
        if any(t == ["N","P","K","CEC","OC","Moisture"])
            y = max(y, realmin);
            yModel = log10(y);
            backTransform = @(z) 10 .^ z;
        else
            yModel = y;
            backTransform = @(z) z;
        end

        cv = cvpartition(height(D), 'KFold', 5);
        yhatPhys = nan(size(y));
        yhatRF = nan(size(y));
        coefMat = [];

        for k = 1:cv.NumTestSets
            tr = training(cv, k);
            te = test(cv, k);

            Xtr = X(tr, :);
            Xte = X(te, :);
            ytr = yModel(tr);

            [B, fitInfo] = lasso(Xtr, ytr, 'CV', 5, 'Alpha', 0.7, 'Standardize', true);
            idx = fitInfo.IndexMinMSE;
            coef = B(:, idx);
            intercept = fitInfo.Intercept(idx);
            predPhys = Xte * coef + intercept;
            yhatPhys(te) = backTransform(predPhys);
            coefMat = [coefMat, [intercept; coef]]; %#ok<AGROW>

            mu = mean(Xtr, 1);
            sig = std(Xtr, 0, 1);
            sig(sig == 0) = 1;
            XtrS = (Xtr - mu) ./ sig;
            XteS = (Xte - mu) ./ sig;

            rf = TreeBagger(300, XtrS, ytr, ...
                'Method', 'regression', ...
                'MinLeafSize', 4, ...
                'NumPredictorsToSample', max(1, round(sqrt(size(XtrS, 2)))));
            yhatRF(te) = backTransform(predict(rf, XteS));
        end

        [physRMSE, physMAE, physR2, physRSS, physTSS] = localRegressionMetrics(y, yhatPhys);
        [rfRMSE, rfMAE, rfR2, rfRSS] = localRegressionMetrics(y, yhatRF);

        coefTbl = table(["Intercept"; featureNames(:)], mean(coefMat, 2), std(coefMat, 0, 2), ...
            'VariableNames', {'Term','MeanCoefficient','StdCoefficient'});
        writetable(coefTbl, fullfile(outDir, sprintf('inverse_coefficients_%s.csv', t)));

        predTbl = table(y, yhatPhys, yhatRF, 'VariableNames', {'Observed','PhysicsInverse','RandomForest'});
        writetable(predTbl, fullfile(outDir, sprintf('inverse_predictions_%s.csv', t)));

        f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
        obs = max(y, realmin);
        pred = max(yhatPhys, realmin);
        scatter(obs, pred, 12, 'filled', 'MarkerFaceAlpha', 0.22);
        hold on;
        lo = min([obs; pred]);
        hi = max([obs; pred]);
        plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.5);
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel(sprintf('Observed %s', t));
        ylabel(sprintf('Predicted %s', t));
        title(sprintf('Inverse model: %s', t));
        txt = sprintf('R^2_{log10} = %.3f | RMSE = %.3g', physR2, physRMSE);
        text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'Margin', 6);
        grid on;
        exportgraphics(f, fullfile(outDir, sprintf('inverse_predicted_vs_actual_%s.png', t)), 'Resolution', 200);
        close(f);

        summaryRows(end+1, 1) = t; %#ok<AGROW>
        summaryN(end+1, 1) = height(D); %#ok<AGROW>
        summaryRMSE(end+1, 1) = physRMSE; %#ok<AGROW>
        summaryMAE(end+1, 1) = physMAE; %#ok<AGROW>
        summaryR2(end+1, 1) = physR2; %#ok<AGROW>

        targetModel = struct();
        targetModel.target = t;
        targetModel.predictors = predictors;
        targetModel.featureNames = featureNames;
        targetModel.coefficients = coefTbl;
        targetModel.metrics = struct( ...
            'Physics_RMSE', physRMSE, 'Physics_MAE', physMAE, 'Physics_R2_log10', physR2, ...
            'Physics_RSS_log10', physRSS, 'Physics_TSS_log10', physTSS, ...
            'RF_RMSE', rfRMSE, 'RF_MAE', rfMAE, 'RF_R2_log10', rfR2, 'RF_RSS_log10', rfRSS, ...
            'N', height(D));
        model.targetModels.(char(t)) = targetModel;
    end

    summaryTbl = table(summaryRows, summaryN, summaryRMSE, summaryMAE, summaryR2, ...
        'VariableNames', {'Target','N','Physics_RMSE','Physics_MAE','Physics_R2_log10'});
    writetable(summaryTbl, fullfile(outDir, 'inverse_model_summary.csv'));

    fid = fopen(fullfile(outDir, 'inverse_model_summary.txt'), 'w');
    fprintf(fid, 'Inverse soil-health model\n');
    fprintf(fid, 'Sensor variables: %s\n', localDisplayList(sensorVars));
    fprintf(fid, 'Context variables: %s\n', localDisplayList(contextVars));
    fprintf(fid, 'Predicted targets: %s\n\n', localDisplayList(summaryRows));
    for i = 1:height(summaryTbl)
        fprintf(fid, '  %s: N=%d, RMSE=%.6g, MAE=%.6g, R^2=%.4f\n', ...
            summaryTbl.Target(i), summaryTbl.N(i), summaryTbl.Physics_RMSE(i), ...
            summaryTbl.Physics_MAE(i), summaryTbl.Physics_R2_log10(i));
    end
    fclose(fid);

    save(fullfile(outDir, 'inverse_soil_health_model.mat'), 'model');
end

function T = localNormalizeVariableNames(T)
    names = string(T.Properties.VariableNames);
    if ismember("H_2O_Vol_", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2O_Vol_;
    end
    if ismember("H_2OVol", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2OVol;
    end
    if ismember("H_2O_Vol", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2O_Vol;
    end
end

function [X, names] = localInverseFeatures(T, sensorVars, contextVars)
    X = [];
    names = strings(0, 1);

    for v = sensorVars
        x = localTransform(v, T.(v));
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "sensor_" + v; %#ok<AGROW>
    end

    for v = contextVars
        x = localTransform(v, T.(v));
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "context_" + v; %#ok<AGROW>
    end

    keyContexts = intersect(["Clay","Silt","Sand","CaCO3"], contextVars, 'stable');
    for sv = sensorVars
        xs = localTransform(sv, T.(sv));
        for cv = keyContexts
            xc = localTransform(cv, T.(cv));
            X = [X, xs .* xc]; %#ok<AGROW>
            names(end+1, 1) = sv + "_x_" + cv; %#ok<AGROW>
        end
    end
end

function x = localTransform(varName, x)
    x = double(x);
    skewed = ["EC","EC_sat","N","P","K","CEC","OC","CaCO3","Coarse","Clay","Moisture"];
    if any(varName == skewed)
        x = log1p(max(x, 0));
    end
end

function [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat)
    keep = isfinite(y) & isfinite(yhat) & y > 0 & yhat > 0;
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

function txt = localDisplayList(x)
    if isempty(x)
        txt = 'none';
    else
        txt = strjoin(cellstr(x), ', ');
    end
end
