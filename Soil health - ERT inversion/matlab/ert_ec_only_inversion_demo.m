clear; clc; close all;

saveArg = true;
useParallel = true;
maxIter = 12;
tolObjective = 1e-3;
tolStep = 1e-3;
targetMisfit = 1e-2;
lineSearchSteps = [1.0, 0.5, 0.25, 0.1, 0.05];

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));
outDir = fullfile(rootDir, 'outputs', 'ert_ec_only_inversion_demo');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Mesh and survey

mesh = make_ert_mesh_fvm_rect( ...
    'XLim', [-30, 30], ...
    'YLim', [-20, 20], ...
    'ZLim', [0, 20], ...
    'Spacing', 10);

survey = make_ert_line_survey( ...
    'XElectrodes', -25:5:25, ...
    'Y', 0, ...
    'Z', 0, ...
    'Current', 1, ...
    'ArrayType', 'wenner');

%% Synthetic true EC from the soil-property forward model

config = make_soil_property_config(mesh);
[mRefProps, refFields] = make_soil_reference_model(config);
[mTrueProps, trueFields] = make_soil_true_model(config, mesh, refFields);

aux = struct();
aux.pH_CaCl2 = 6.4 * ones(config.nCells, 1);
aux.pH_H2O = 7.0 * ones(config.nCells, 1);
aux.CaCO3 = 2.0 * ones(config.nCells, 1);
aux.Coarse = 5.0 * ones(config.nCells, 1);

ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

[sigmaTrue, ~, ~] = predict_soil_ec_physics(mTrueProps, config, aux, ecModel);
Ktrue = assemble_ert_stiffness(mesh, sigmaTrue, ...
    'UseParallel', useParallel, ...
    'ElementType', mesh.elementType, ...
    'ShowProgress', true);
resultTrue = solve_ert_forward(mesh, survey, Ktrue);
dObs = resultTrue.voltage;

%% Initial EC model

sigmaStart = predict_soil_ec_physics(mRefProps, config, aux, ecModel);
sigmaStart = sigmaStart(:);
uRef = log10(max(sigmaStart, realmin));
uCurrent = uRef;
sigmaCurrent = sigmaStart;

dataStd = max(1e-3, 0.03 * max(abs(dObs), 1e-3));
Wd = spdiags(1 ./ dataStd, 0, numel(dObs), numel(dObs));
ecWeights = default_ec_inversion_weights();

Kcurrent = assemble_ert_stiffness(mesh, sigmaCurrent, ...
    'UseParallel', useParallel, ...
    'ElementType', mesh.elementType, ...
    'ShowProgress', true);
resultCurrent = solve_ert_forward(mesh, survey, Kcurrent);
dPred = resultCurrent.voltage;
dStart = dPred;
objectiveCurrent = evaluate_ec_inversion_objective( ...
    uCurrent, dObs, dPred, Wd, mesh, uRef, ecWeights);

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
Jert = [];

%% Iterative EC-only inversion

for iter = 1:maxIter
    fprintf('EC-only iteration %d of %d\n', iter, maxIter);

    Jert = build_ert_jacobian_fd(mesh, survey, sigmaCurrent, resultCurrent, ...
        'UseParallel', useParallel, ...
        'ShowProgress', true, ...
        'RelativePerturbation', 0.02);
    Ju = Jert * spdiags(log(10) * sigmaCurrent, 0, numel(sigmaCurrent), numel(sigmaCurrent));

    [R, rhsReg] = build_ec_regularization(mesh, uCurrent, uRef, ecWeights);
    rData = dObs - dPred;
    [deltaU, ~] = solve_property_update(Ju, rData, Wd, R, rhsReg);

    accepted = false;
    bestTrial = struct();
    for stepLength = lineSearchSteps
        uTrial = uCurrent + stepLength * deltaU;
        sigmaTrial = 10 .^ uTrial;
        Ktrial = assemble_ert_stiffness(mesh, sigmaTrial, ...
            'UseParallel', useParallel, ...
            'ElementType', mesh.elementType, ...
            'ShowProgress', true);
        resultTrial = solve_ert_forward(mesh, survey, Ktrial);
        dTrial = resultTrial.voltage;
        objectiveTrial = evaluate_ec_inversion_objective( ...
            uTrial, dObs, dTrial, Wd, mesh, uRef, ecWeights);

        if isfinite(objectiveTrial.total) && objectiveTrial.total < objectiveCurrent.total
            accepted = true;
            bestTrial.u = uTrial;
            bestTrial.sigma = sigmaTrial;
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
    iterationLog(iter).stepNorm = norm(deltaU) / max(norm(uCurrent), eps);
    iterationLog(iter).stepLength = lastStepLength;
    iterationLog(iter).accepted = accepted;

    if ~accepted
        exitReason = "lineSearchFailed";
        lastDeltaU = deltaU;
        lastStepLength = 0;
        break;
    end

    relativeObjectiveDrop = (objectiveCurrent.total - bestTrial.objective.total) / max(objectiveCurrent.total, eps);
    relativeStepNorm = norm(bestTrial.stepLength * deltaU) / max(norm(uCurrent), eps);

    uCurrent = bestTrial.u;
    sigmaCurrent = bestTrial.sigma;
    resultCurrent = bestTrial.result;
    dPred = bestTrial.dPred;
    objectiveCurrent = bestTrial.objective;
    lastDeltaU = deltaU;
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
summary.objectiveBefore = evaluate_ec_inversion_objective(uRef, dObs, dStart, Wd, mesh, uRef, ecWeights).total;
summary.objectiveAfter = objectiveCurrent.total;
summary.stepNorm = norm(lastStepLength * lastDeltaU) / max(norm(uCurrent), eps);
summary.lastStepLength = lastStepLength;
summary.iterationsCompleted = numel(iterationLog);
summary.exitReason = exitReason;
summary.normalizedDataMisfitAfter = objectiveCurrent.normalizedDataMisfit;

fig = plot_ec_only_inversion_demo(mesh, sigmaTrue, sigmaStart, sigmaUpdated, dObs, dPred, iterationLog, targetMisfit, summary);
figPseudo = plot_ert_pseudosection(survey, dObs, dPred);

if saveArg
    save(fullfile(outDir, 'ert_ec_only_inversion_demo.mat'), ...
        'mesh', 'survey', 'config', 'aux', 'ecModel', 'ecWeights', ...
        'mRefProps', 'mTrueProps', 'trueFields', 'refFields', ...
        'sigmaTrue', 'sigmaStart', 'sigmaUpdated', 'uRef', 'uCurrent', ...
        'dObs', 'dUpdated', 'Jert', 'iterationLog', 'summary', 'exitReason');
    exportgraphics(fig, fullfile(outDir, 'ert_ec_only_inversion_demo.png'), 'Resolution', 220);
    exportgraphics(figPseudo, fullfile(outDir, 'ert_ec_only_pseudosection.png'), 'Resolution', 220);
end
