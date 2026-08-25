clear; clc;

% Estimate local soil-property influence on the direct Random Forest EC model.
% Because the RF model is nonparametric, sensitivities are finite-difference
% perturbation responses around the latest staged property inversion result.

relativePerturbation = 0.01; % 1 percent increase
absoluteSteps = struct( ...
    'N', 1.0, ...
    'P', 1.0, ...
    'K', 10.0, ...
    'OC', 0.1, ...
    'Clay', 1.0, ...
    'Sand', 1.0, ...
    'Silt', 1.0, ...
    'pH_CaCl2', 0.1, ...
    'pH_H2O', 0.1, ...
    'CaCO3', 1.0, ...
    'Coarse', 1.0);

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(scriptDir, 'scr')));
addpath(fullfile(rootDir, 'soil_health_ec_models'));

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
propertyFile = fullfile(outDir, 'staged_property_inversion_output.mat');
rfModelFile = fullfile(rootDir, 'soil_health_ec_models', ...
    'direct_ec_random_forest_direct_ec_rf_model.mat');
if ~isfile(propertyFile)
    error('Missing property inversion output: %s', propertyFile);
end
if ~isfile(rfModelFile)
    error('Missing random forest model: %s', rfModelFile);
end

S = load(propertyFile);
rf = load(rfModelFile);
model = rf.model;

config = S.config;
if isfield(S, 'fieldsProp')
    fieldsBase = S.fieldsProp;
elseif isfield(S, 'fieldsStart')
    fieldsBase = S.fieldsStart;
else
    error('Property output does not contain fieldsProp or fieldsStart.');
end

if isfield(S, 'activeMask') && numel(S.activeMask) == config.nCells
    activeMask = logical(S.activeMask(:));
else
    activeMask = true(config.nCells, 1);
end

predictorNames = string(model.predictorNames);
baseTable = localBuildPredictorTable(fieldsBase, config.nCells, predictorNames);
ecBase = localPredictDirectEcRandomForest(baseTable, model);

rows = cell(numel(predictorNames), 1);
for i = 1:numel(predictorNames)
    name = predictorNames(i);
    xBase = baseTable.(char(name));
    absStep = localGetStep(absoluteSteps, name);

    relTable = baseTable;
    relTable.(char(name)) = max(xBase .* (1 + relativePerturbation), 0);
    ecRel = localPredictDirectEcRandomForest(relTable, model);
    relResp = localEcFractionChange(ecRel, ecBase, activeMask);
    actualRelFrac = median((relTable.(char(name)) - xBase) ./ max(abs(xBase), eps), 'omitnan');
    elasticity = relResp ./ max(actualRelFrac, eps);

    absTable = baseTable;
    absTable.(char(name)) = max(xBase + absStep, 0);
    ecAbs = localPredictDirectEcRandomForest(absTable, model);
    absResp = localEcFractionChange(ecAbs, ecBase, activeMask);

    rows{i} = {char(name), median(xBase(activeMask), 'omitnan'), ...
        100 * median(relResp, 'omitnan'), 100 * mean(relResp, 'omitnan'), ...
        100 * max(abs(relResp), [], 'omitnan'), ...
        median(elasticity, 'omitnan'), mean(elasticity, 'omitnan'), ...
        absStep, 100 * median(absResp, 'omitnan'), ...
        100 * mean(absResp, 'omitnan'), 100 * max(abs(absResp), [], 'omitnan')};
end

T = cell2table(vertcat(rows{:}), 'VariableNames', { ...
    'Property', 'MedianValue', ...
    'MedianPctEcChangeFrom1PctPropertyIncrease', ...
    'MeanPctEcChangeFrom1PctPropertyIncrease', ...
    'MaxAbsPctEcChangeFrom1PctPropertyIncrease', ...
    'MedianElasticityPctEcPerPctProperty', ...
    'MeanElasticityPctEcPerPctProperty', ...
    'AbsoluteStep', ...
    'MedianPctEcChangeFromAbsoluteStep', ...
    'MeanPctEcChangeFromAbsoluteStep', ...
    'MaxAbsPctEcChangeFromAbsoluteStep'});

if isfield(model, 'importance') && istable(model.importance)
    importance = model.importance;
    importance.Predictor = string(importance.Predictor);
    T.Property = string(T.Property);
    T = outerjoin(T, importance, 'LeftKeys', 'Property', 'RightKeys', 'Predictor', ...
        'MergeKeys', false, 'Type', 'left');
    T.Predictor = [];
    T = movevars(T, 'OOBImportance', 'After', 'Property');
end

csvOut = fullfile(outDir, 'random_forest_ec_parameter_sensitivity.csv');
txtOut = fullfile(outDir, 'random_forest_ec_parameter_sensitivity.txt');
writetable(T, csvOut);

fid = fopen(txtOut, 'w');
cleanup = onCleanup(@() fclose(fid));
localPrintTable(fid, T, relativePerturbation, propertyFile, rfModelFile);
localPrintTable(1, T, relativePerturbation, propertyFile, rfModelFile);

fprintf('\nSaved random forest sensitivity table:\n  %s\n  %s\n', csvOut, txtOut);

function T = localBuildPredictorTable(fields, nCells, predictorNames)
    T = table();
    for name = predictorNames
        key = char(name);
        switch key
            case {'N', 'P', 'K', 'OC', 'Clay', 'Sand', 'Silt'}
                if ~isfield(fields, key)
                    error('Missing field required by RF model: %s', key);
                end
                T.(key) = double(fields.(key)(:));
            case 'pH_CaCl2'
                T.(key) = 6.4 * ones(nCells, 1);
            case 'pH_H2O'
                T.(key) = 7.0 * ones(nCells, 1);
            case 'CaCO3'
                T.(key) = 2.0 * ones(nCells, 1);
            case 'Coarse'
                T.(key) = 5.0 * ones(nCells, 1);
            otherwise
                error('Unsupported RF predictor: %s', key);
        end
    end
end

function step = localGetStep(steps, name)
    key = char(name);
    if isfield(steps, key)
        step = steps.(key);
    else
        step = 1.0;
    end
end

function response = localEcFractionChange(ecPert, ecBase, activeMask)
    response = (ecPert(activeMask) - ecBase(activeMask)) ./ max(abs(ecBase(activeMask)), realmin);
end

function ecHat = localPredictDirectEcRandomForest(newData, model)
    predictorNames = string(model.predictorNames);
    T = newData(:, predictorNames);
    for v = string(model.skewedPredictors)
        if ismember(v, predictorNames)
            T.(char(v)) = log1p(max(T.(char(v)), 0));
        end
    end
    Xraw = T{:, cellstr(predictorNames)};
    X = (Xraw - model.mu) ./ model.sigma;
    ecHat = exp(predict(model.rf, X)) - 1;
end

function localPrintTable(fid, T, relativePerturbation, propertyFile, rfModelFile)
    fprintf(fid, 'Random Forest soil -> EC local sensitivity\n');
    fprintf(fid, '  Property source: %s\n', propertyFile);
    fprintf(fid, '  RF model: %s\n', rfModelFile);
    fprintf(fid, '  Relative perturbation: %.3g%% property increase\n\n', 100 * relativePerturbation);
    fprintf(fid, '  Elasticity is approximately (%% EC change) / (%% property change).\n');
    fprintf(fid, '  Positive means increasing the property raises EC; negative lowers EC.\n');
    fprintf(fid, '  OOBImportance is the forest training-time permutation importance.\n\n');

    hasImportance = any(strcmp(T.Properties.VariableNames, 'OOBImportance'));
    if hasImportance
        fprintf(fid, '%-11s %12s %12s %14s %14s %14s\n', ...
            'property', 'median', 'OOBImp', 'medElastic', 'absStep', 'medAbsStep%');
        fprintf(fid, '%s\n', repmat('-', 1, 85));
        for i = 1:height(T)
            fprintf(fid, '%-11s %12.4g %12.4g %14.4g %14.4g %14.4g\n', ...
                string(T.Property(i)), T.MedianValue(i), T.OOBImportance(i), ...
                T.MedianElasticityPctEcPerPctProperty(i), T.AbsoluteStep(i), ...
                T.MedianPctEcChangeFromAbsoluteStep(i));
        end
    else
        fprintf(fid, '%-11s %12s %14s %14s %14s\n', ...
            'property', 'median', 'medElastic', 'absStep', 'medAbsStep%');
        fprintf(fid, '%s\n', repmat('-', 1, 72));
        for i = 1:height(T)
            fprintf(fid, '%-11s %12.4g %14.4g %14.4g %14.4g\n', ...
                string(T.Property(i)), T.MedianValue(i), ...
                T.MedianElasticityPctEcPerPctProperty(i), T.AbsoluteStep(i), ...
                T.MedianPctEcChangeFromAbsoluteStep(i));
        end
    end
end
