clear; clc; close all;

saveArg = true;
useParallel = true;
maxIter = 8;
tolObjective = 1e-4;
tolStep = 1e-3;
targetMisfit = 1e-1;
lineSearchSteps = [1.0, 0.5, 0.25, 0.1, 0.05];
refinementFactor = 5;
% EC_uniform_high | EC_uniform_low | 
% EC_gradient_horizontal | EC_gradient_vertical
startModelType = 'EC_uniform_high'; 
parameterizationType = "background_plus_anomaly"; % background_plus_anomaly | cellwise

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));
outDir = fullfile(rootDir, 'outputs', 'ert_ec_only_inversion_demo_2d');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
startModelTag = regexprep(char(lower(startModelType)), '[^a-z0-9_]+', '_');

%% Mesh and survey
% Invert on a coarser x-z mesh, but solve the electric potentials on a
% finer forward mesh.
electrodeX = -25:5:25;
xLim = [-120, 120];
yLim = [0, 5];
zLimCoarse = [0, 58];
zLimFine = [0, 58];

% Finer around the electrode line and near surface, coarser away from the
% main sensing region.
xCoarse = make_piecewise_axis(xLim, [-40, 40], [20, 10, 20]);
xFine = make_piecewise_axis(xLim, [-40, 40], [10, 5, 10]);
zCoarse = make_piecewise_axis(zLimCoarse, [6, 18], [2, 4, 8]);
zFine = make_piecewise_axis(zLimFine, [6, 18], [1, 2, 4]);

mesh = make_ert_mesh_fvm_rect( ...
    'XLim', xLim, ...
    'YLim', yLim, ...
    'ZLim', zLimCoarse, ...
    'XCoords', xCoarse, ...
    'YCoords', yLim, ...
    'ZCoords', zCoarse);

forwardMesh = make_ert_mesh_fvm_rect( ...
    'XLim', xLim, ...
    'YLim', yLim, ...
    'ZLim', zLimFine, ...
    'XCoords', xFine, ...
    'YCoords', yLim, ...
    'ZCoords', zFine);

[fineFromCoarse, coarseFromFine] = make_structured_cell_mapping(mesh, forwardMesh);

survey = make_ert_line_survey( ...
    'XElectrodes', electrodeX, ...
    'Y', 0, ...
    'Z', 0, ...
    'Current', 1, ...
    'ArrayType', ["wenner", "dipole-dipole", "schlumberger"]);

%% Synthetic true EC from the soil-property forward model

config = make_soil_property_config(mesh);
configForward = make_soil_property_config(forwardMesh);
[mRefProps, refFields] = make_soil_reference_model(config);
[mTruePropsForward, trueFieldsForward] = make_soil_true_model(configForward, forwardMesh, localReferenceFields(configForward));

aux = struct();
aux.pH_CaCl2 = 6.4 * ones(configForward.nCells, 1);
aux.pH_H2O = 7.0 * ones(configForward.nCells, 1);
aux.CaCO3 = 2.0 * ones(configForward.nCells, 1);
aux.Coarse = 5.0 * ones(configForward.nCells, 1);

auxCoarse = struct();
auxCoarse.pH_CaCl2 = coarseFromFine * aux.pH_CaCl2;
auxCoarse.pH_H2O = coarseFromFine * aux.pH_H2O;
auxCoarse.CaCO3 = coarseFromFine * aux.CaCO3;
auxCoarse.Coarse = coarseFromFine * aux.Coarse;

ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

[sigmaTrueForward, ~, ~] = predict_soil_ec_physics(mTruePropsForward, configForward, aux, ecModel);
sigmaTrue = coarseFromFine * sigmaTrueForward;
trueFields = localAverageFieldsToCoarse(trueFieldsForward, coarseFromFine, config.names);

Ktrue = assemble_ert_stiffness(forwardMesh, sigmaTrueForward, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true);
resultTrue = solve_ert_forward(forwardMesh, survey, Ktrue);
dObs = resultTrue.voltage;

%% Initial EC model

sigmaStart = localBuildInitialEcModel(startModelType, mesh, forwardMesh, sigmaTrueForward, coarseFromFine);
sigmaStart = sigmaStart(:);
uStart = log10(max(sigmaStart, realmin));

if parameterizationType == "background_plus_anomaly"
    param = make_background_anomaly_parameterization(mesh);
    [uBgRef, uAnomRef] = decompose_u_background_anomaly(mesh, uStart);
    uBgCurrent = uBgRef;
    uAnomCurrent = uAnomRef;
    pCurrent = [uBgCurrent; uAnomCurrent];
    uRef = param.G * [uBgRef; uAnomRef];
    uCurrent = param.G * pCurrent;
    ecWeights = default_ec_background_anomaly_weights();
else
    param = [];
    uRef = uStart;
    uCurrent = uRef;
    ecWeights = default_ec_inversion_weights();
end
sigmaCurrent = 10 .^ uCurrent;

dataStd = max(1e-3, 0.03 * max(abs(dObs), 1e-3));
Wd = spdiags(1 ./ dataStd, 0, numel(dObs), numel(dObs));

sigmaForwardCurrent = 10 .^ (fineFromCoarse * uCurrent);
Kcurrent = assemble_ert_stiffness(forwardMesh, sigmaForwardCurrent, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true);
resultCurrent = solve_ert_forward(forwardMesh, survey, Kcurrent);
dPred = resultCurrent.voltage;
dStart = dPred;
if parameterizationType == "background_plus_anomaly"
    objectiveCurrent = evaluate_ec_background_anomaly_objective( ...
        dObs, dPred, Wd, mesh, param, uBgCurrent, uAnomCurrent, uBgRef, uAnomRef, ecWeights);
else
    objectiveCurrent = evaluate_ec_inversion_objective( ...
        uCurrent, dObs, dPred, Wd, mesh, uRef, ecWeights);
end

iterationLog = repmat(struct( ...
    'iter', 0, ...
    'objective', 0, ...
    'phiData', 0, ...
    'phiReg', 0, ...
    'normalizedDataMisfit', 0, ...
    'stepNorm', 0, ...
    'stepLength', 0, ...
    'accepted', false), maxIter, 1);
exitReason = "maxIter";
lastDeltaU = zeros(size(uCurrent));
lastStepLength = 0;
Ju = [];

%% Iterative EC-only inversion

for iter = 1:maxIter
    fprintf('EC-only 2D iteration %d of %d\n', iter, maxIter);

    Ju = build_ert_jacobian_fd_u(forwardMesh, survey, uCurrent, fineFromCoarse, resultCurrent, ...
        'UseParallel', useParallel, ...
        'ShowProgress', true, ...
        'RelativePerturbation', 0.02, ...
        'AbsolutePerturbation', 1e-3);

    if parameterizationType == "background_plus_anomaly"
        Jsolve = Ju * param.G;
        [R, rhsReg] = build_ec_background_anomaly_regularization( ...
            mesh, param, uBgCurrent, uAnomCurrent, uBgRef, uAnomRef, ecWeights);
    else
        Jsolve = Ju;
        [R, rhsReg] = build_ec_regularization(mesh, uCurrent, uRef, ecWeights);
    end
    rData = dObs - dPred;
    [deltaSolve, ~] = solve_property_update(Jsolve, rData, Wd, R, rhsReg);

    accepted = false;
    bestTrial = struct();
    for stepLength = lineSearchSteps
        if parameterizationType == "background_plus_anomaly"
            pTrial = pCurrent + stepLength * deltaSolve;
            uBgTrial = pTrial(param.idxBg);
            uAnomTrial = pTrial(param.idxAnom);
            uTrial = param.G * pTrial;
        else
            uTrial = uCurrent + stepLength * deltaSolve;
        end
        sigmaTrial = 10 .^ uTrial;
        sigmaTrialForward = 10 .^ (fineFromCoarse * uTrial);
        Ktrial = assemble_ert_stiffness(forwardMesh, sigmaTrialForward, ...
            'UseParallel', useParallel, ...
            'ElementType', forwardMesh.elementType, ...
            'ShowProgress', true);
        resultTrial = solve_ert_forward(forwardMesh, survey, Ktrial);
        dTrial = resultTrial.voltage;
        if parameterizationType == "background_plus_anomaly"
            objectiveTrial = evaluate_ec_background_anomaly_objective( ...
                dObs, dTrial, Wd, mesh, param, uBgTrial, uAnomTrial, uBgRef, uAnomRef, ecWeights);
        else
            objectiveTrial = evaluate_ec_inversion_objective( ...
                uTrial, dObs, dTrial, Wd, mesh, uRef, ecWeights);
        end

        if isfinite(objectiveTrial.total) && objectiveTrial.total < objectiveCurrent.total
            accepted = true;
            if parameterizationType == "background_plus_anomaly"
                bestTrial.p = pTrial;
                bestTrial.uBg = uBgTrial;
                bestTrial.uAnom = uAnomTrial;
            end
            bestTrial.u = uTrial;
            bestTrial.sigma = sigmaTrial;
            bestTrial.sigmaForward = sigmaTrialForward;
            bestTrial.result = resultTrial;
            bestTrial.dPred = dTrial;
            bestTrial.objective = objectiveTrial;
            bestTrial.stepLength = stepLength;
            break;
        end
    end

    iterationLog(iter).iter = iter;
    iterationLog(iter).objective = objectiveCurrent.total;
    iterationLog(iter).phiData = objectiveCurrent.phiData;
    iterationLog(iter).phiReg = objectiveCurrent.phiReg;
    iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
    iterationLog(iter).stepNorm = norm(deltaSolve) / max(norm(uCurrent), eps);
    iterationLog(iter).stepLength = lastStepLength;
    iterationLog(iter).accepted = accepted;

    if ~accepted
        exitReason = "lineSearchFailed";
        lastDeltaU = deltaSolve;
        lastStepLength = 0;
        break;
    end

    relativeObjectiveDrop = (objectiveCurrent.total - bestTrial.objective.total) / max(objectiveCurrent.total, eps);
    relativeStepNorm = norm(bestTrial.stepLength * deltaSolve) / max(norm(uCurrent), eps);

    if parameterizationType == "background_plus_anomaly"
        pCurrent = bestTrial.p;
        uBgCurrent = bestTrial.uBg;
        uAnomCurrent = bestTrial.uAnom;
    end
    uCurrent = bestTrial.u;
    sigmaCurrent = bestTrial.sigma;
    sigmaForwardCurrent = bestTrial.sigmaForward;
    resultCurrent = bestTrial.result;
    dPred = bestTrial.dPred;
    objectiveCurrent = bestTrial.objective;
    lastDeltaU = deltaSolve;
    lastStepLength = bestTrial.stepLength;

    iterationLog(iter).objective = objectiveCurrent.total;
    iterationLog(iter).phiData = objectiveCurrent.phiData;
    iterationLog(iter).phiReg = objectiveCurrent.phiReg;
    iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
    iterationLog(iter).stepNorm = relativeStepNorm;
    iterationLog(iter).stepLength = bestTrial.stepLength;

    if objectiveCurrent.normalizedDataMisfit <= targetMisfit
        exitReason = "targetMisfit";
        break;
    end
    if relativeObjectiveDrop < tolObjective
        exitReason = "smallObjectiveImprovement";
        break;
    end
    if relativeStepNorm < tolStep
        exitReason = "smallStep";
        break;
    end
end

iterationLog = iterationLog([iterationLog.iter] > 0);
sigmaUpdated = sigmaCurrent;
dUpdated = dPred;

summary = struct();
summary.nCells = numel(sigmaTrue);
summary.nMeasurements = numel(dObs);
summary.dataMisfitBefore = norm(Wd * (dObs - dStart));
summary.dataMisfitAfter = norm(Wd * (dObs - dUpdated));
if parameterizationType == "background_plus_anomaly"
    summary.objectiveBefore = evaluate_ec_background_anomaly_objective( ...
        dObs, dStart, Wd, mesh, param, uBgRef, uAnomRef, uBgRef, uAnomRef, ecWeights).total;
else
    summary.objectiveBefore = evaluate_ec_inversion_objective(uRef, dObs, dStart, Wd, mesh, uRef, ecWeights).total;
end
summary.objectiveAfter = objectiveCurrent.total;
summary.stepNorm = norm(lastStepLength * lastDeltaU) / max(norm(uCurrent), eps);
summary.lastStepLength = lastStepLength;
summary.iterationsCompleted = numel(iterationLog);
summary.exitReason = exitReason;
summary.normalizedDataMisfitAfter = objectiveCurrent.normalizedDataMisfit;
summary.refinementFactor = refinementFactor;
summary.parameterizationType = parameterizationType;
summary.startModelType = startModelType;


%% Plotting

fig = plot_ec_only_inversion_demo( ...
    mesh, sigmaTrue, sigmaStart, sigmaUpdated, dObs, dPred, iterationLog, targetMisfit, summary, ...
    'ForwardMesh', forwardMesh, ...
    'SigmaTrueForward', sigmaTrueForward, ...
    'Survey', survey, ...
    'TolObjective', tolObjective);
% [figPotential, selectedSource] = plot_selected_source_potential( ...
%     forwardMesh, survey, resultTrue, resultCurrent, 'SignedLogMinAbs', 1e-3);
% figPseudo = plot_ert_pseudosection(survey, dObs, dPred);

if saveArg
    % save(fullfile(outDir, ['ert_ec_only_inversion_demo_2d_' startModelTag '.mat']), ...
    %     'mesh', 'forwardMesh', 'fineFromCoarse', 'coarseFromFine', ...
    %     'survey', 'config', 'configForward', 'aux', 'auxCoarse', 'ecModel', 'ecWeights', ...
    %     'mRefProps', 'mTruePropsForward', 'trueFields', 'trueFieldsForward', 'refFields', ...
    %     'sigmaTrue', 'sigmaTrueForward', 'sigmaStart', 'sigmaUpdated', 'sigmaForwardCurrent', ...
    %     'uRef', 'uCurrent', 'dObs', 'dUpdated', 'Ju', 'iterationLog', 'summary', ...
    %     'exitReason', 'resultTrue', 'resultCurrent', 'selectedSource');
    exportgraphics(fig, fullfile(outDir, ['ert_ec_only_inversion_demo_2d_' startModelTag '.png']), 'Resolution', 220);
    % exportgraphics(figPotential, fullfile(outDir, ['ert_selected_source_potential_2d_' startModelTag '.png']), 'Resolution', 220);
    % exportgraphics(figPseudo, fullfile(outDir, 'ert_ec_only_pseudosection_2d.png'), 'Resolution', 220);
end


%% Functions

function refFields = localReferenceFields(configLocal)
%LOCALREFERENCEFIELDS Homogeneous reference used to seed the synthetic truth.

    refFields = struct();
    refFields.N = 30 * ones(configLocal.nCells, 1);
    refFields.P = 20 * ones(configLocal.nCells, 1);
    refFields.K = 120 * ones(configLocal.nCells, 1);
    refFields.OC = 2.5 * ones(configLocal.nCells, 1);
    refFields.Moisture = 18 * ones(configLocal.nCells, 1);
    refFields.Clay = 30 * ones(configLocal.nCells, 1);
    refFields.Sand = 45 * ones(configLocal.nCells, 1);
    refFields.Silt = 25 * ones(configLocal.nCells, 1);
end

function fieldsCoarse = localAverageFieldsToCoarse(fieldsFine, coarseFromFine, names)
%LOCALAVERAGEFIELDSTOCOARSE Average fine synthetic fields back to the inversion grid.

    fieldsCoarse = struct();
    for i = 1:numel(names)
        name = char(names(i));
        fieldsCoarse.(name) = coarseFromFine * fieldsFine.(name);
    end
end

function sigmaStart = localBuildInitialEcModel(startModelType, coarseMesh, forwardMesh, sigmaTrueForward, coarseFromFine)
%LOCALBUILDINITIALECMODEL Build one of several deterministic EC starting models.

    sigmaTrueCoarse = coarseFromFine * sigmaTrueForward(:);
    [highVal, lowVal] = localEstimateBackgroundLevels(forwardMesh, sigmaTrueForward(:));
    centers = coarseMesh.elementCenters;

    switch string(startModelType)
        case "EC_uniform_high"
            sigmaStart = highVal * ones(coarseMesh.nElements, 1);

        case "EC_uniform_low"
            sigmaStart = lowVal * ones(coarseMesh.nElements, 1);

        case "EC_gradient_horizontal"
            x = centers(:, 1);
            xNorm = (x - min(x)) ./ max(max(x) - min(x), eps);
            sigmaStart = lowVal + (highVal - lowVal) * xNorm;

        case "EC_gradient_vertical"
            z = centers(:, 3);
            zNorm = (z - min(z)) ./ max(max(z) - min(z), eps);
            sigmaStart = lowVal + (highVal - lowVal) * zNorm;

        otherwise
            error('Unsupported startModelType: %s', startModelType);
    end

    sigmaStart = max(min(sigmaStart, max(sigmaTrueCoarse) * 5), min(sigmaTrueCoarse) * 0.1);
end

function [highVal, lowVal] = localEstimateBackgroundLevels(forwardMesh, sigmaTrueForward)
%LOCALESTIMATEBACKGROUNDLEVELS Estimate shallow-high and layered-low background EC.

    c = forwardMesh.elementCenters;
    x = c(:, 1);
    z = c(:, 3);

    backgroundMask = abs(x) >= 0.35 * max(abs(x));
    if ~any(backgroundMask)
        backgroundMask = true(size(x));
    end

    [zLevels, ~, zGroup] = unique(round(z, 10), 'stable');
    zMeans = nan(size(zLevels));
    for i = 1:numel(zLevels)
        mask = backgroundMask & (zGroup == i);
        if any(mask)
            zMeans(i) = mean(sigmaTrueForward(mask));
        end
    end

    valid = isfinite(zMeans);
    if ~any(valid)
        highVal = mean(sigmaTrueForward);
        lowVal = mean(sigmaTrueForward);
        return;
    end

    firstValid = find(valid, 1, 'first');
    highVal = zMeans(firstValid);
    lowVal = min(zMeans(valid));
end
