function model = fit_physics_ec_extract_model(T, options)
% Fit a physics-informed EC model for aqueous soil extracts.
%
% Model:
%   sigma_pred = beta0 + a*sigma_sol + b*clayFrac*h(I,pH) + c*OC + d*CEC + e*CaCO3
%
% where:
%   h(I,pH) = log1p(I/I0)^gamma * exp(eta*(pH - pHref))
%
% The fit is performed against log10(EC) residuals so that relative errors
% matter more than absolute errors at high EC.
%
% Required target column:
%   EC_sat  measured EC in S/m
%
% Required ion columns for best use:
%   Ca_mgL, Mg_mgL, Na_mgL, K_mgL, Cl_mgL, SO4_mgL, NO3_mgL, HCO3_mgL
%
% If those are unavailable, the model automatically falls back to a
% proxy-chemistry mode using pH, K, P, N, CEC, CaCO3, and OC.

    arguments
        T table
        options.TargetName string = "EC_sat"
        options.TemperatureC double = 25
        options.I0_molL double = 1e-3
        options.pHref double = 7
        options.SaveModelFile string = "physics_ec_extract_model.mat"
    end

    if ~ismember(options.TargetName, string(T.Properties.VariableNames))
        error('Target column "%s" was not found.', options.TargetName);
    end
    paths = codex_paths('physics_ec_extract');
    outDir = paths.analysisDir;

    features = compute_extract_chemistry_features(T, 'TemperatureC', options.TemperatureC);
    y = T.(options.TargetName);

    D = [features, table(y, 'VariableNames', {'EC_measured_Sm'})];
    keep = all(isfinite(D{:,:}), 2) & D.EC_measured_Sm > 0;
    D = D(keep, :);

    p0 = [ ...
        0.001, ...  % beta0 baseline S/m
        1.0,   ...  % a solution-conductivity multiplier
        0.01,  ...  % b clay/surface multiplier
        0.6,   ...  % gamma ionic-strength nonlinearity
        0.05,  ...  % eta pH sensitivity
        0.0,   ...  % c OC term
        0.0,   ...  % d CEC term
        0.0];       % e CaCO3 term

    mdl = fitnlm(D, @(p, X) localPhysicsPredict(p, X, options.I0_molL, options.pHref), ...
        log10(D.EC_measured_Sm), p0);

    sigma_pred = 10 .^ predict(mdl, D);
    [rmse, mae, r2, rss, tss] = localRegressionMetrics(D.EC_measured_Sm, sigma_pred);

    model = struct();
    model.type = "physics_extract_ec";
    model.fit = mdl;
    model.options = options;
    model.featureNames = string(features.Properties.VariableNames);
    model.chemistryMode = ternary(any(D.chemistryMode > 0), "ion-based", "proxy-based");
    model.metrics = struct('N', height(D), 'RMSE', rmse, 'MAE', mae, ...
        'R2_log10', r2, 'RSS_log10', rss, 'TSS_log10', tss);

    saveFile = char(options.SaveModelFile);
    if ~contains(saveFile, filesep)
        saveFile = fullfile(outDir, saveFile);
    end
    save(saveFile, 'model');

    coeff = mdl.Coefficients;
    writetable(coeff, fullfile(outDir, 'physics_ec_extract_coefficients.csv'));

    fid = fopen(fullfile(outDir, 'physics_ec_extract_summary.txt'), 'w');
    fprintf(fid, 'Physics-informed extract EC model\n');
    fprintf(fid, 'Chemistry mode: %s\n', model.chemistryMode);
    fprintf(fid, 'Rows used: %d\n', model.metrics.N);
    fprintf(fid, 'RMSE = %.6g S/m\n', model.metrics.RMSE);
    fprintf(fid, 'MAE = %.6g S/m\n', model.metrics.MAE);
    fprintf(fid, 'RSS in log10 space = %.6g\n', model.metrics.RSS_log10);
    fprintf(fid, 'TSS in log10 space = %.6g\n', model.metrics.TSS_log10);
    fprintf(fid, 'R^2 in log10 space = %.4f\n', model.metrics.R2_log10);
    fclose(fid);
end

function ylog = localPhysicsPredict(p, X, I0, pHref)
    sigmaSol = X.sigma_sol_Sm;
    I = max(X.ionicStrength_molL, 0);
    pH = X.pH;
    clayFrac = max(X.clayFrac, 0);
    oc = max(X.OC_gkg, 0);
    cec = max(X.CEC, 0);
    caco3 = max(X.CaCO3, 0);

    h = (log1p(I ./ I0) .^ p(4)) .* exp(p(5) .* (pH - pHref));
    sigmaPred = p(1) + p(2).*sigmaSol + p(3).*clayFrac.*h + p(6).*oc + p(7).*cec + p(8).*caco3;
    sigmaPred = max(sigmaPred, realmin);
    ylog = log10(sigmaPred);
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

function [rmse, mae, r2, rss, tss] = localRegressionMetrics(y, yhat)
    y = y(:);
    yhat = yhat(:);
    rmse = sqrt(mean((y - yhat).^2));
    mae = mean(abs(y - yhat));
    yLog = log10(max(y, realmin));
    yhatLog = log10(max(yhat, realmin));
    tss = sum((yLog - mean(yLog)).^2);
    rss = sum((yLog - yhatLog).^2);
    r2 = 1 - rss / tss;
end
