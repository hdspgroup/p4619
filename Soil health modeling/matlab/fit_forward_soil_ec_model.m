function model = fit_forward_soil_ec_model(T, options)
% Fit a forward soil EC model:
%   state variables + context variables -> EC
%
% State variables (targets you eventually care about):
%   N, P, K, CEC, OC, Moisture
%
% Context variables (soil modifiers / conditioning variables):
%   Clay, Silt, Sand, Coarse, pH_CaCl2, pH_H2O, CaCO3
%
% This function supports missing optional variables. If Moisture or CEC are
% absent they are omitted from the active fit, but moisture is still tracked
% as an expected placeholder for future datasets.
%
% The fit is done in log10(EC) space using a hybrid feature design:
%   - direct state terms
%   - direct context terms
%   - selected state x context interactions
%   - Random Forest benchmark for comparison
%
% Current physics-informed model form used in this script:
%
%   Let y = log10(EC)
%
%   y_hat = beta0
%         + sum_i beta_state_i * g_i(state_i)
%         + sum_j beta_context_j * h_j(context_j)
%         + sum_{i,j in selected pairs} beta_ij * g_i(state_i) * h_j(context_j)
%
% where:
%   state_i   in {N, P, K, CEC, OC, Moisture if available}
%   context_j in {Clay, Silt, Sand, Coarse, pH_CaCl2, pH_H2O, CaCO3}
%
% transforms:
%   g_i(x) = log1p(x) for skewed positive state variables
%   h_j(x) = log1p(x) for skewed positive context variables
%   h_j(x) = x otherwise
%
% selected interaction pairs currently include each state variable with:
%   {Clay, Silt, Sand, pH_H2O, pH_CaCl2, CaCO3}
%
% final EC prediction:
%   EC_hat = 10^(y_hat)
%
% Example:
%   T = readtable('my_forward_training_data.csv');
%   model = fit_forward_soil_ec_model(T, 'TargetName', 'EC_sat');

    arguments
        T table
        options.TargetName string = "EC_sat"
        options.AnalysisName string = "forward_soil_ec"
        options.TargetUnits string = "S/m"
        options.ConvertECdSmToSm logical = false
    end

    if ~ismember(options.TargetName, string(T.Properties.VariableNames))
        error('Target column "%s" was not found.', options.TargetName);
    end

    paths = codex_paths(options.AnalysisName);
    outDir = paths.analysisDir;

    expectedStateVars = ["N","P","K","CEC","OC","Moisture"];
    moistureAliases = ["Moisture","H_2O_Vol_","H_2OVol","H_2O_Vol","H_2OVol","WaterContent","theta_v"];
    stateCandidates = ["N","P","K","CEC","OC", moistureAliases];
    contextCandidates = ["Clay","Silt","Sand","Coarse","pH_CaCl2","pH_H2O","CaCO3"];

    stateVars = intersect(stateCandidates, string(T.Properties.VariableNames), 'stable');
    contextVars = intersect(contextCandidates, string(T.Properties.VariableNames), 'stable');

    for a = moistureAliases
        if ismember(a, string(T.Properties.VariableNames)) && a ~= "Moisture"
            T.Moisture = T.(a);
            stateVars(stateVars == a) = "Moisture";
        end
    end
    stateVars = unique(stateVars, 'stable');
    missingExpectedState = setdiff(expectedStateVars, stateVars, 'stable');

    if options.ConvertECdSmToSm
        T.(options.TargetName) = T.(options.TargetName) .* 0.1;
        options.TargetUnits = "S/m (from dS/m input)";
    end

    useVars = [stateVars, contextVars, options.TargetName];
    D = T(:, useVars);
    D = D(all(~ismissing(D), 2), :);
    D = D(isfinite(D.(options.TargetName)) & D.(options.TargetName) > 0, :);

    [Xphys, physNames, meta] = localForwardFeatures(D, stateVars, contextVars);
    y = D.(options.TargetName);
    yLog = log10(max(y, realmin));
    modelEquation = localForwardModelEquation();

    cv = cvpartition(height(D), 'KFold', 5);
    yhatPhys = nan(size(y));
    yhatRF = nan(size(y));
    coefMat = [];

    for k = 1:cv.NumTestSets
        tr = training(cv, k);
        te = test(cv, k);

        Xtr = Xphys(tr, :);
        Xte = Xphys(te, :);
        ytr = yLog(tr);

        [B, fitInfo] = lasso(Xtr, ytr, 'CV', 5, 'Alpha', 0.7, 'Standardize', true);
        idx = fitInfo.IndexMinMSE;
        coef = B(:, idx);
        intercept = fitInfo.Intercept(idx);
        predPhys = Xte * coef + intercept;
        yhatPhys(te) = 10 .^ predPhys;
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
        yhatRF(te) = 10 .^ predict(rf, XteS);
    end

    [physRMSE, physMAE, physR2, physRSS, tss] = localRegressionMetrics(y, yhatPhys);
    [rfRMSE, rfMAE, rfR2, rfRSS] = localRegressionMetrics(y, yhatRF);

    summaryTbl = table( ...
        ["Physics-informed forward"; "Random Forest benchmark"], ...
        [height(D); height(D)], ...
        [physRMSE; rfRMSE], ...
        [physMAE; rfMAE], ...
        [physRSS; rfRSS], ...
        [physR2; rfR2], ...
        'VariableNames', {'Model','N','RMSE','MAE','RSS_log10','R2_log10'});
    writetable(summaryTbl, fullfile(outDir, 'forward_model_summary.csv'));

    coefNames = ["Intercept"; physNames(:)];
    coefTbl = table(coefNames, mean(coefMat, 2), std(coefMat, 0, 2), ...
        'VariableNames', {'Term','MeanCoefficient','StdCoefficient'});
    writetable(coefTbl, fullfile(outDir, 'forward_model_coefficients.csv'));

    predTbl = table(y, yhatPhys, yhatRF, 'VariableNames', {'Observed','PhysicsForward','RandomForest'});
    writetable(predTbl, fullfile(outDir, 'forward_model_predictions.csv'));

    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
    scatter(y, yhatPhys, 12, 'filled', 'MarkerFaceAlpha', 0.22);
    hold on;
    lo = min([y; yhatPhys]);
    hi = max([y; yhatPhys]);
    plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel(sprintf('Observed %s', options.TargetName));
    ylabel(sprintf('Predicted %s', options.TargetName));
    title('Forward EC model: physics-informed fit');
    txt = sprintf('R^2_{log10} = %.3f | RMSE = %.3g', physR2, physRMSE);
    text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 6, 'FontSize', 14, 'FontWeight', 'bold');
    set(gca,'FontSize',14);
    grid on;
    exportgraphics(f, fullfile(outDir, 'forward_model_predicted_vs_actual.png'), 'Resolution', 200);
    close(f);

    fid = fopen(fullfile(outDir, 'forward_model_summary.txt'), 'w');
    fprintf(fid, 'Forward soil EC model\n');
    fprintf(fid, 'Target: %s\n', options.TargetName);
    fprintf(fid, 'Units: %s\n', options.TargetUnits);
    fprintf(fid, 'Rows used: %d\n', height(D));
    fprintf(fid, 'State variables: %s\n', strjoin(cellstr(stateVars), ', '));
    fprintf(fid, 'Missing optional state placeholders: %s\n', localDisplayList(missingExpectedState));
    fprintf(fid, 'Context variables: %s\n', strjoin(cellstr(contextVars), ', '));
    fprintf(fid, 'TSS in log10 space: %.6g\n\n', tss);
    fprintf(fid, 'Physics-informed forward: RMSE = %.6g, MAE = %.6g, RSS = %.6g, R^2 = %.4f\n', ...
        physRMSE, physMAE, physRSS, physR2);
    fprintf(fid, 'Random Forest benchmark: RMSE = %.6g, MAE = %.6g, RSS = %.6g, R^2 = %.4f\n', ...
        rfRMSE, rfMAE, rfRSS, rfR2);
    fprintf(fid, '\nInterpretation:\n');
    fprintf(fid, '  State variables are the soil-health quantities we care about.\n');
    fprintf(fid, '  Context variables condition how those state variables translate into EC.\n');
    fprintf(fid, '  Clay belongs in the context block because it changes the EC response.\n');
    fclose(fid);

    fid = fopen(fullfile(outDir, 'forward_model_equation.txt'), 'w');
    fprintf(fid, '%s', modelEquation);
    fclose(fid);

    model = struct();
    model.type = "forward_soil_ec";
    model.expectedStateVars = expectedStateVars;
    model.stateVars = stateVars;
    model.missingStateVars = missingExpectedState;
    model.contextVars = contextVars;
    model.featureNames = physNames;
    model.featureMeta = meta;
    model.equation = modelEquation;
    model.metrics = struct( ...
        'Physics_RMSE', physRMSE, 'Physics_MAE', physMAE, 'Physics_R2_log10', physR2, ...
        'RF_RMSE', rfRMSE, 'RF_MAE', rfMAE, 'RF_R2_log10', rfR2, ...
        'TSS_log10', tss, 'N', height(D));
    model.coefficients = coefTbl;

    save(fullfile(outDir, 'forward_soil_ec_model.mat'), 'model');
end

function [X, names, meta] = localForwardFeatures(T, stateVars, contextVars)
    n = height(T);
    X = [];
    names = strings(0, 1);

    for v = stateVars
        x = T.(v);
        x = localTransform(v, x);
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "state_" + v; %#ok<AGROW>
    end

    for v = contextVars
        x = T.(v);
        x = localTransform(v, x);
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "context_" + v; %#ok<AGROW>
    end

    % Add selected interactions: every state variable with the main structural modifiers.
    keyContexts = intersect(["Clay","Silt","Sand","pH_H2O","pH_CaCl2","CaCO3"], contextVars, 'stable');
    for sv = stateVars
        xs = localTransform(sv, T.(sv));
        for cv = keyContexts
            xc = localTransform(cv, T.(cv));
            X = [X, xs .* xc]; %#ok<AGROW>
            names(end+1, 1) = sv + "_x_" + cv; %#ok<AGROW>
        end
    end

    meta = struct();
    meta.stateVars = stateVars;
    meta.contextVars = contextVars;
    meta.keyContexts = keyContexts;
    meta.nRows = n;
end

function x = localTransform(varName, x)
    x = double(x);
    skewed = ["N","P","K","CEC","OC","CaCO3","Coarse","Clay","Moisture"];
    if any(varName == skewed)
        x = log1p(max(x, 0));
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

function txt = localDisplayList(x)
    if isempty(x)
        txt = 'none';
    else
        txt = strjoin(cellstr(x), ', ');
    end
end

function txt = localForwardModelEquation()
    lines = strings(0,1);
    lines(end+1) = "Current physics-informed forward soil EC model";
    lines(end+1) = "";
    lines(end+1) = "Purpose";
    lines(end+1) = "Predict EC from soil-health state variables plus soil context variables.";
    lines(end+1) = "";
    lines(end+1) = "Response";
    lines(end+1) = "y = log10(EC)";
    lines(end+1) = "";
    lines(end+1) = "Model form";
    lines(end+1) = "y_hat = beta0";
    lines(end+1) = "      + sum_i beta_state_i * g_i(state_i)";
    lines(end+1) = "      + sum_j beta_context_j * h_j(context_j)";
    lines(end+1) = "      + sum_(i,j in selected pairs) beta_ij * g_i(state_i) * h_j(context_j)";
    lines(end+1) = "";
    lines(end+1) = "State variables";
    lines(end+1) = "N, P, K, CEC, OC, Moisture if available";
    lines(end+1) = "";
    lines(end+1) = "Context variables";
    lines(end+1) = "Clay, Silt, Sand, Coarse, pH_CaCl2, pH_H2O, CaCO3";
    lines(end+1) = "";
    lines(end+1) = "Transforms";
    lines(end+1) = "For skewed positive variables:";
    lines(end+1) = "  g(x) = log1p(x)";
    lines(end+1) = "  h(x) = log1p(x)";
    lines(end+1) = "For non-skewed variables:";
    lines(end+1) = "  use x directly";
    lines(end+1) = "";
    lines(end+1) = "Selected interactions";
    lines(end+1) = "Each state variable is multiplied by the key context variables:";
    lines(end+1) = "  Clay, Silt, Sand, pH_H2O, pH_CaCl2, CaCO3";
    lines(end+1) = "";
    lines(end+1) = "Prediction";
    lines(end+1) = "EC_hat = 10^(y_hat)";
    lines(end+1) = "";
    lines(end+1) = "Interpretation";
    lines(end+1) = "State variables represent soil-health quantities of interest.";
    lines(end+1) = "Context variables condition how those state variables translate into conductivity.";
    lines(end+1) = "Clay belongs in the context block because it modifies the EC response rather than acting only as a health target.";
    txt = strjoin(lines, newline);
end
