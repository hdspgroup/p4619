clear; clc;

% Estimate local soil-property influence on the physics-informed EC model.
% The reported elasticity is approximately:
%     (% change in EC) / (% change in property)
% evaluated around the latest staged soil-property inversion result.

relativePerturbation = 0.01; % 1 percent increase
absoluteSteps = struct( ...
    'N', 1.0, ...
    'P', 1.0, ...
    'K', 10.0, ...
    'OC', 0.1, ...
    'Moisture', 1.0, ...
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

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
propertyFile = fullfile(outDir, 'staged_property_inversion_output.mat');
if ~isfile(propertyFile)
    error('Missing property inversion output: %s', propertyFile);
end

S = load(propertyFile);
ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

mesh = S.mesh;
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

auxBase = localMakeAux(config.nCells);
mBase = pack_soil_fields(fieldsBase, config);
sigmaBase = predict_soil_ec_physics(mBase, config, auxBase, ecModel);

variableNames = unique([string(config.names), "pH_CaCl2", "pH_H2O", "CaCO3", "Coarse"], 'stable');
rows = cell(numel(variableNames), 1);
for i = 1:numel(variableNames)
    name = variableNames(i);
    xBase = localGetVariable(fieldsBase, auxBase, name);
    absStep = localGetStep(absoluteSteps, name);

    [sigmaRel, actualRelFrac] = localPerturbAndPredict( ...
        fieldsBase, auxBase, config, ecModel, name, relativePerturbation, "relative");
    relResp = localEcFractionChange(sigmaRel, sigmaBase, activeMask);
    elasticity = relResp ./ max(actualRelFrac, eps);

    [sigmaAbs, actualAbsStep] = localPerturbAndPredict( ...
        fieldsBase, auxBase, config, ecModel, name, absStep, "absolute");
    absResp = localEcFractionChange(sigmaAbs, sigmaBase, activeMask);

    rows{i} = {char(name), median(xBase(activeMask), 'omitnan'), ...
        100 * median(relResp, 'omitnan'), 100 * mean(relResp, 'omitnan'), ...
        100 * max(abs(relResp), [], 'omitnan'), ...
        median(elasticity, 'omitnan'), mean(elasticity, 'omitnan'), ...
        actualAbsStep, 100 * median(absResp, 'omitnan'), ...
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

csvOut = fullfile(outDir, 'soil_ec_parameter_sensitivity.csv');
txtOut = fullfile(outDir, 'soil_ec_parameter_sensitivity.txt');
writetable(T, csvOut);

fid = fopen(txtOut, 'w');
cleanup = onCleanup(@() fclose(fid));
localPrintTable(fid, T, relativePerturbation, propertyFile);
localPrintTable(1, T, relativePerturbation, propertyFile);

fprintf('\nSaved sensitivity table:\n  %s\n  %s\n', csvOut, txtOut);

function aux = localMakeAux(nCells)
    aux = struct();
    aux.pH_CaCl2 = 6.4 * ones(nCells, 1);
    aux.pH_H2O = 7.0 * ones(nCells, 1);
    aux.CaCO3 = 2.0 * ones(nCells, 1);
    aux.Coarse = 5.0 * ones(nCells, 1);
end

function x = localGetVariable(fields, aux, name)
    key = char(name);
    if isfield(fields, key)
        x = fields.(key);
    elseif isfield(aux, key)
        x = aux.(key);
    else
        error('Unknown variable: %s', key);
    end
    x = double(x(:));
end

function step = localGetStep(steps, name)
    key = char(name);
    if isfield(steps, key)
        step = steps.(key);
    else
        step = 1.0;
    end
end

function [sigmaPert, actualPerturbation] = localPerturbAndPredict(fieldsBase, auxBase, config, ecModel, name, amount, mode)
    fieldsPert = fieldsBase;
    auxPert = auxBase;
    key = char(name);
    x = localGetVariable(fieldsBase, auxBase, name);
    switch lower(string(mode))
        case "relative"
            xPert = x .* (1 + amount);
            actualPerturbation = median((xPert - x) ./ max(abs(x), eps), 'omitnan');
        case "absolute"
            xPert = x + amount;
            actualPerturbation = amount;
        otherwise
            error('Unsupported perturbation mode: %s', mode);
    end

    if isfield(fieldsPert, key)
        fieldsPert.(key) = xPert;
    elseif isfield(auxPert, key)
        auxPert.(key) = xPert;
    else
        error('Unknown variable: %s', key);
    end

    mPert = pack_soil_fields(fieldsPert, config);
    sigmaPert = predict_soil_ec_physics(mPert, config, auxPert, ecModel);
end

function response = localEcFractionChange(sigmaPert, sigmaBase, activeMask)
    response = (sigmaPert(activeMask) - sigmaBase(activeMask)) ./ max(abs(sigmaBase(activeMask)), realmin);
end

function localPrintTable(fid, T, relativePerturbation, propertyFile)
    fprintf(fid, 'Physics-informed soil -> EC local sensitivity\n');
    fprintf(fid, '  Source: %s\n', propertyFile);
    fprintf(fid, '  Relative perturbation: %.3g%% property increase\n\n', 100 * relativePerturbation);
    fprintf(fid, '  Elasticity is approximately (%% EC change) / (%% property change).\n');
    fprintf(fid, '  Positive means increasing the property raises EC; negative lowers EC.\n\n');
    fprintf(fid, '%-11s %12s %14s %14s %14s %14s\n', ...
        'property', 'median', 'medElastic', 'meanElastic', 'absStep', 'medAbsStep%');
    fprintf(fid, '%s\n', repmat('-', 1, 83));
    for i = 1:height(T)
        fprintf(fid, '%-11s %12.4g %14.4g %14.4g %14.4g %14.4g\n', ...
            string(T.Property{i}), T.MedianValue(i), ...
            T.MedianElasticityPctEcPerPctProperty(i), ...
            T.MeanElasticityPctEcPerPctProperty(i), ...
            T.AbsoluteStep(i), ...
            T.MedianPctEcChangeFromAbsoluteStep(i));
    end
end
