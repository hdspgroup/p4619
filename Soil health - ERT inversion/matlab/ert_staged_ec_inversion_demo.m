clear; clc; close all;

saveArg = true;
useParallel = true;
progressMode = "text";
forwardMeshType = "fvm_octree_style"; % fvm_octree_style | fvm_rect_refined
lateralBuffer = 15;

% Stage plan controls
nGeometricGroups = 3;
backgroundSchedule = {1, 1:2};
maxIterBackground = 3;
maxIterAnomaly = 5;
lineSearchBackground = 1*[2.0, 1.0, 0.5, 0.25];
lineSearchAnomaly = 2*[2.0, 1.0, 0.5, 0.25, 0.1];
plotSliceZ = 3.0;

% Horizontal e-folding radius of the N anomaly; the model uses radius^2 internally.
truthSettings = struct( ...
    'layerCenterZ', 2, ...
    'layerWidthDenom', 0.4, ...
    'anomalyCenterZ', 0.5, ...
    'anomalyXYRadius', 5.0, ...
    'anomalyZDenom', 0.28, ...
    'nBackground', 30, ...
    'nAnomalyAmplitude', 6, ...
    'nAnomalyFloor', 29, ...
    'clayBackground', 5, ...
    'clayLayerAmplitude', 0, ...
    'clayFloor', 5, ...
    'sandBackground', 70, ...
    'sandLayerAmplitude', 0, ...
    'sandFloor', 5);
truthSettings.anomalyXYDenom = truthSettings.anomalyXYRadius.^2;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

priorCsvPath = fullfile(rootDir, 'data', 'prior_samples', 'soil_property_prior_samples_aligned_dense.csv');
weightsCsvPath = fullfile(rootDir, 'matlab', 'property_regularization_weights.csv');
surveyDir = fullfile(rootDir, 'data', 'electrode_arrays');
surveyFile = localFindMostRecentSurveyFile(surveyDir);
surveyData = load(surveyFile, 'survey', 'summary');
survey = assign_survey_geometric_factor_groups(surveyData.survey, nGeometricGroups);

fprintf('\nStaged inversion survey grouping\n');
for g = 1:nGeometricGroups
    mask = survey.grouping.groupIndex == g;
    fprintf('  Group %d: %d measurements | %s\n', g, nnz(mask), survey.grouping.groupLabels(g));
end
fprintf('\n');

%% Meshes
surveyCenter = mean(survey.electrodes(:, 1:2), 1);
surveyRadius = max(vecnorm(survey.electrodes(:, 1:2) - surveyCenter, 2, 2));
outerRadius = surveyRadius + lateralBuffer;

sb = max(25, surveyRadius);
orad = outerRadius;
xfc = [5, 10, 20, sb, orad];
xff = sort([2.5, xfc, mean([xfc(1:end-1); xfc(2:end)], 1)]);
zfc = [0, 1, 2, 4, 6];
zff = sort([zfc, mean([zfc(1:end-1); zfc(2:end)], 1)]);

xCoarse = surveyCenter(1) + [-fliplr(xfc), 0, xfc];
yCoarse = surveyCenter(2) + [-fliplr(xfc), 0, xfc];
zCoarse = zfc;
xFine = surveyCenter(1) + [-fliplr(xff), 0, xff];
yFine = surveyCenter(2) + [-fliplr(xff), 0, xff];
zFine = zff;

mesh = make_ert_mesh_fvm_rect( ...
    'XLim', [xCoarse(1), xCoarse(end)], ...
    'YLim', [yCoarse(1), yCoarse(end)], ...
    'ZLim', [zCoarse(1), zCoarse(end)], ...
    'XCoords', xCoarse, ...
    'YCoords', yCoarse, ...
    'ZCoords', zCoarse);

switch lower(forwardMeshType)
    case "fvm_rect_refined"
        forwardMesh = make_ert_mesh_fvm_rect( ...
            'XLim', [xFine(1), xFine(end)], ...
            'YLim', [yFine(1), yFine(end)], ...
            'ZLim', [zFine(1), zFine(end)], ...
            'XCoords', xFine, ...
            'YCoords', yFine, ...
            'ZCoords', zFine);
    case "fvm_octree_style"
        forwardMesh = make_ert_mesh_fvm_octree_style( ...
            'Electrodes', survey.electrodes, ...
            'CenterXY', surveyCenter, ...
            'Radius', outerRadius, ...
            'CoreHalfWidth', max(10, 0.7 * surveyRadius), ...
            'MinCellSizeXY', 2.5, ...
            'MaxCellSizeXY', 10.0, ...
            'SurfaceRefineDepth', 3.0, ...
            'MinCellSizeZ', 0.5, ...
            'MaxCellSizeZ', 2.0, ...
            'ZMax', zFine(end), ...
            'OuterPadding', lateralBuffer);
    otherwise
        error('Unsupported forwardMeshType: %s', forwardMeshType);
end

[fineFromCoarse, coarseFromFine] = make_structured_cell_mapping(mesh, forwardMesh);
localPrintProblemSummary(mesh, forwardMesh, survey, surveyFile, useParallel, forwardMeshType);

%% Soil truth and priors
ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

config = make_soil_property_config(mesh);
configForward = make_soil_property_config(forwardMesh);
[~, refFields] = make_soil_reference_model(config);
[~, refFieldsForward] = make_soil_reference_model(configForward);
[mTrue, fieldsTrue] = make_soil_true_model(config, mesh, refFields, truthSettings);
[mTrueForward, fieldsTrueForward] = make_soil_true_model(configForward, forwardMesh, refFieldsForward, truthSettings);

if contains(string(priorCsvPath), "soil_property_prior_samples_aligned_dense.csv", 'IgnoreCase', true)
    localRefreshAlignedPriorCsv(priorCsvPath, refFields, ecModel, truthSettings);
end
aux = localMakeAux(config.nCells);
auxForward = localMakeAux(configForward.nCells);

[sigmaTrueForward, ~] = predict_soil_ec_physics(mTrueForward, configForward, auxForward, ecModel);
sigmaTrue = coarseFromFine * sigmaTrueForward;

Ktrue = assemble_ert_stiffness(forwardMesh, sigmaTrueForward, ...
    'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, 'ProgressMode', progressMode);
resultTrue = solve_ert_forward(forwardMesh, survey, Ktrue);
dObs = resultTrue.voltage;

%% Stage 1: background-only EC inversion
bgWeights = default_ec_background_anomaly_weights(weightsCsvPath);
localPrintParsedEcWeights(bgWeights, weightsCsvPath);
param = make_background_anomaly_parameterization(mesh, 'BackgroundMode', bgWeights.backgroundMode);
fprintf('  EC background parameterization: %s (%d background coefficient(s))\n', ...
    bgWeights.backgroundMode, param.nBg);
sigmaStart = mean(sigmaTrue) * ones(mesh.nElements, 1);
uStart = log10(max(sigmaStart, realmin));
[uBgRef, ~] = decompose_u_background_anomaly(mesh, uStart, param);
uBgCurrent = uBgRef;
uAnomZero = zeros(mesh.nElements, 1);

iterationLogBackground = struct([]);
for s = 1:numel(backgroundSchedule)
    groupsHere = backgroundSchedule{s};
    mask = ismember(survey.grouping.groupIndex, groupsHere);
    surveyStage = subset_ert_survey(survey, mask);
    fprintf('Background stage %d using groups [%s] with %d measurements\n', ...
        s, sprintf('%d ', groupsHere), surveyStage.nMeasurements);

    [uBgCurrent, stageLog] = localRunBackgroundStage( ...
        mesh, forwardMesh, surveyStage, fineFromCoarse, dObs(mask), ...
        survey, dObs, ...
        uBgCurrent, uBgRef, param, bgWeights, maxIterBackground, ...
        lineSearchBackground, useParallel, progressMode);
    iterationLogBackground = [iterationLogBackground; stageLog(:)]; %#ok<AGROW>
end

uStage1 = param.B * uBgCurrent;
sigmaStage1 = 10 .^ uStage1;
sigmaStage1Forward = 10 .^ (fineFromCoarse * uStage1);
Kstage1 = assemble_ert_stiffness(forwardMesh, sigmaStage1Forward, ...
    'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, 'ProgressMode', progressMode);
resultStage1 = solve_ert_forward(forwardMesh, survey, Kstage1);

%% Stage 2: anomaly EC inversion around the recovered background
uBgStage2 = uBgCurrent;
uAnomStage2 = zeros(mesh.nElements, 1);

[uBgStage2, uAnomStage2, sigmaStage2, resultStage2, iterationLogEc, ecSummary] = localRunAnomalyStage( ...
    mesh, forwardMesh, survey, fineFromCoarse, dObs, param, ...
    uBgStage2, uAnomStage2, uBgCurrent, zeros(mesh.nElements,1), ...
    bgWeights, maxIterAnomaly, lineSearchAnomaly, useParallel, progressMode);
localPrintAnomalyWeightGuidance(ecSummary, bgWeights, fullfile(outDir, 'stage2_anomaly_weight_guidance.txt'));

uStage2 = param.G * [uBgStage2; uAnomStage2];
sigmaStage2 = sigmaStage2(:);

%% Plot and save
ecSummary.trueModelMeta = localStageTrueMeta(forwardMesh, truthSettings, plotSliceZ);
ecSummary.cellSensitivity = [];
ecSummary.parameterizationType = "background_plus_anomaly_staged";
ecSummary.forwardMeshType = forwardMeshType;
ecSummary.surveyFile = surveyFile;
ecSummary.surveySummary = surveyData.summary;

ecTrueSlices = localBuildAnalyticEcSlices(forwardMesh, refFieldsForward, truthSettings, ecModel, auxForward, plotSliceZ);
figEc = plot_staged_ec_inversion_demo( ...
    mesh, survey, sigmaStage1, sigmaStage2, ecTrueSlices, ...
    dObs, resultStage1.voltage, resultStage2.voltage, ...
    iterationLogBackground, iterationLogEc, ecSummary);

if saveArg
    exportgraphics(figEc, fullfile(outDir, 'ert_staged_ec_inversion_demo.png'), 'Resolution', 220);
    save(fullfile(outDir, 'staged_ec_inversion_output.mat'), ...
        'mesh', 'forwardMesh', 'survey', 'surveyFile', 'truthSettings', ...
        'sigmaTrue', 'sigmaTrueForward', 'sigmaStage1', 'sigmaStage2', ...
        'fieldsTrue', 'fieldsTrueForward', 'mTrue', 'mTrueForward', ...
        'dObs', 'iterationLogBackground', 'iterationLogEc', ...
        'resultTrue', 'resultStage1', 'resultStage2', ...
        'ecSummary', 'plotSliceZ', 'forwardMeshType', 'weightsCsvPath', 'priorCsvPath');
    save(fullfile(outDir, 'ert_staged_ec_property_inversion_demo.mat'), ...
        'mesh', 'forwardMesh', 'survey', 'surveyFile', 'truthSettings', ...
        'sigmaTrue', 'sigmaTrueForward', 'sigmaStage1', 'sigmaStage2', ...
        'fieldsTrue', 'fieldsTrueForward', 'mTrue', 'mTrueForward', ...
        'dObs', 'iterationLogBackground', 'iterationLogEc', ...
        'resultTrue', 'resultStage1', 'resultStage2', ...
        'ecSummary', 'plotSliceZ', 'forwardMeshType', 'weightsCsvPath', 'priorCsvPath');
end

%% Local functions

function [uBgCurrent, stageLog] = localRunBackgroundStage(mesh, forwardMesh, surveyStage, fineFromCoarse, dObsStage, ...
        surveyFull, dObsFull, uBgCurrent, uBgRef, param, weights, maxIter, lineSearchSteps, useParallel, progressMode)

    nData = numel(dObsStage);
    dataStd = max(1e-3, 0.03 * max(abs(dObsStage), 1e-3));
    Wd = spdiags(1 ./ dataStd, 0, nData, nData);
    stageLog = repmat(struct('iter',0,'objective',0,'phiData',0,'phiReg',0,'normalizedDataMisfit',0,'stepNorm',0,'stepLength',0,'accepted',false), maxIter, 1);

    uAnomZero = zeros(mesh.nElements, 1);
    uCurrent = param.B * uBgCurrent;
    sigmaForwardCurrent = 10 .^ (fineFromCoarse * uCurrent);
    Kcurrent = assemble_ert_stiffness(forwardMesh, sigmaForwardCurrent, ...
        'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
        'ShowProgress', true, 'ProgressMode', progressMode);
    resultCurrent = solve_ert_forward(forwardMesh, surveyStage, Kcurrent);
    dPred = resultCurrent.voltage;
    fullDataStd = max(1e-3, 0.03 * max(abs(dObsFull), 1e-3));
    WdFull = spdiags(1 ./ fullDataStd, 0, numel(dObsFull), numel(dObsFull));
    resultCurrentFull = solve_ert_forward(forwardMesh, surveyFull, Kcurrent);
    dPredFull = resultCurrentFull.voltage;

    for iter = 1:maxIter
        Ju = build_ert_jacobian_fd_u(forwardMesh, surveyStage, uCurrent, fineFromCoarse, resultCurrent, ...
            'UseParallel', useParallel, 'ShowProgress', true, 'ProgressMode', progressMode, ...
            'RelativePerturbation', 0.02, 'AbsolutePerturbation', 1e-3);
        Jbg = Ju * param.B;
        [Rbg, rhsBg] = localBuildBackgroundRegularization(param, uBgCurrent, uBgRef, weights);
        rData = dObsStage - dPred;
        [deltaBg, ~] = solve_property_update(Jbg, rData, Wd, Rbg, rhsBg);

        phiDataCurrent = sum((Wd * rData).^2);
        phiDataCurrentFull = sum((WdFull * (dObsFull - dPredFull)).^2);
        phiRegCurrent = localBackgroundRegularizationValue(param, uBgCurrent, uBgRef, weights);
        objectiveCurrent = phiDataCurrent + phiRegCurrent;

        accepted = false;
        best = struct();
        for stepLength = lineSearchSteps
            uBgTrial = uBgCurrent + stepLength * deltaBg;
            uTrial = param.B * uBgTrial;
            sigmaTrialForward = 10 .^ (fineFromCoarse * uTrial);
            Ktrial = assemble_ert_stiffness(forwardMesh, sigmaTrialForward, ...
                'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
                'ShowProgress', true, 'ProgressMode', progressMode);
            resultTrial = solve_ert_forward(forwardMesh, surveyStage, Ktrial);
            dTrial = resultTrial.voltage;
            rTrial = dObsStage - dTrial;
            phiDataTrial = sum((Wd * rTrial).^2);
            resultTrialFull = solve_ert_forward(forwardMesh, surveyFull, Ktrial);
            dTrialFull = resultTrialFull.voltage;
            phiDataTrialFull = sum((WdFull * (dObsFull - dTrialFull)).^2);
            phiRegTrial = localBackgroundRegularizationValue(param, uBgTrial, uBgRef, weights);
            objectiveTrial = phiDataTrial + phiRegTrial;
            if objectiveTrial < objectiveCurrent
                accepted = true;
                best.uBg = uBgTrial;
                best.u = uTrial;
                best.dPred = dTrial;
                best.result = resultTrial;
                best.dPredFull = dTrialFull;
                best.objective = objectiveTrial;
                best.phiData = phiDataTrial;
                best.phiDataFull = phiDataTrialFull;
                best.phiReg = phiRegTrial;
                best.stepLength = stepLength;
                break;
            end
        end

        stageLog(iter).iter = iter;
        stageLog(iter).accepted = accepted;
        stageLog(iter).objective = phiDataCurrent + phiRegCurrent;
        stageLog(iter).phiData = phiDataCurrentFull;
        stageLog(iter).phiReg = phiRegCurrent;
        stageLog(iter).normalizedDataMisfit = phiDataCurrentFull / max(1, numel(dObsFull));
        stageLog(iter).stepNorm = norm(deltaBg) / max(norm(uBgCurrent), eps);

        if ~accepted
            break;
        end

        uBgCurrent = best.uBg;
        uCurrent = best.u;
        dPred = best.dPred;
        dPredFull = best.dPredFull;
        resultCurrent = best.result;
        stageLog(iter).objective = best.objective;
        stageLog(iter).phiData = best.phiDataFull;
        stageLog(iter).phiReg = best.phiReg;
        stageLog(iter).normalizedDataMisfit = best.phiDataFull / max(1, numel(dObsFull));
        stageLog(iter).stepLength = best.stepLength;
    end

    stageLog = stageLog([stageLog.iter] > 0);
end

function [uBgCurrent, uAnomCurrent, sigmaCurrent, resultCurrent, iterationLog, summary] = localRunAnomalyStage( ...
        mesh, forwardMesh, survey, fineFromCoarse, dObs, param, ...
        uBgCurrent, uAnomCurrent, uBgRef, uAnomRef, weights, maxIter, lineSearchSteps, useParallel, progressMode)

    dataStd = max(1e-3, 0.03 * max(abs(dObs), 1e-3));
    Wd = spdiags(1 ./ dataStd, 0, numel(dObs), numel(dObs));
    uBgFixed = uBgCurrent;
    uCurrent = param.B * uBgFixed + uAnomCurrent;
    sigmaCurrent = 10 .^ uCurrent;
    sigmaForwardCurrent = 10 .^ (fineFromCoarse * uCurrent);
    Kcurrent = assemble_ert_stiffness(forwardMesh, sigmaForwardCurrent, ...
        'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
        'ShowProgress', true, 'ProgressMode', progressMode);
    resultCurrent = solve_ert_forward(forwardMesh, survey, Kcurrent);
    dPred = resultCurrent.voltage;
    objectiveCurrent = localEvaluateAnomalyOnlyObjective(dObs, dPred, Wd, mesh, param, uAnomCurrent, uAnomRef, weights);
    objectiveStart = objectiveCurrent;
    rawResidualStart = dObs - dPred;

    iterationLog = repmat(struct('iter',0,'objective',0,'phiData',0,'phiReg',0,'normalizedDataMisfit',0,'stepNorm',0,'stepLength',0,'accepted',false), maxIter, 1);
    exitReason = "maxIter";
    lastStepLength = 0;
    lastDelta = zeros(size(uAnomCurrent));
    JuFinal = [];
    firstDeltaGuidance = struct();

    for iter = 1:maxIter
        Ju = build_ert_jacobian_fd_u(forwardMesh, survey, uCurrent, fineFromCoarse, resultCurrent, ...
            'UseParallel', useParallel, 'ShowProgress', true, 'ProgressMode', progressMode, ...
            'RelativePerturbation', 0.02, 'AbsolutePerturbation', 1e-3);
        JuFinal = Ju;
        Jsolve = Ju;
        [R, rhsReg] = localBuildAnomalyOnlyRegularization(mesh, param, uAnomCurrent, uAnomRef, weights);
        rData = dObs - dPred;
        [deltaAnom, ~] = solve_property_update(Jsolve, rData, Wd, R, rhsReg);
        if iter == 1
            firstDeltaGuidance = localBuildAnomalyWeightGuidance(mesh, param, deltaAnom, uAnomRef, objectiveStart.phiData);
        end

        accepted = false;
        best = struct();
        for stepLength = lineSearchSteps
            uAnomTrial = uAnomCurrent + stepLength * deltaAnom;
            uTrial = param.B * uBgFixed + uAnomTrial;
            sigmaTrial = 10 .^ uTrial;
            sigmaTrialForward = 10 .^ (fineFromCoarse * uTrial);
            Ktrial = assemble_ert_stiffness(forwardMesh, sigmaTrialForward, ...
                'UseParallel', useParallel, 'ElementType', forwardMesh.elementType, ...
                'ShowProgress', true, 'ProgressMode', progressMode);
            resultTrial = solve_ert_forward(forwardMesh, survey, Ktrial);
            dTrial = resultTrial.voltage;
            objectiveTrial = localEvaluateAnomalyOnlyObjective(dObs, dTrial, Wd, mesh, param, uAnomTrial, uAnomRef, weights);
            if objectiveTrial.total < objectiveCurrent.total
                accepted = true;
                best.uAnom = uAnomTrial;
                best.u = uTrial;
                best.sigma = sigmaTrial;
                best.result = resultTrial;
                best.dPred = dTrial;
                best.objective = objectiveTrial;
                best.stepLength = stepLength;
                break;
            end
        end

        iterationLog(iter).iter = iter;
        iterationLog(iter).objective = objectiveCurrent.total;
        iterationLog(iter).phiData = objectiveCurrent.phiData;
        iterationLog(iter).phiReg = objectiveCurrent.phiReg;
        iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
        iterationLog(iter).stepNorm = norm(deltaAnom) / max(norm(uAnomCurrent), 1);
        iterationLog(iter).accepted = accepted;

        if ~accepted
            exitReason = "lineSearchFailed";
            lastDelta = deltaAnom;
            break;
        end

        uAnomCurrent = best.uAnom;
        uCurrent = best.u;
        sigmaCurrent = best.sigma;
        resultCurrent = best.result;
        dPred = best.dPred;
        objectiveCurrent = best.objective;
        lastDelta = deltaAnom;
        lastStepLength = best.stepLength;

        iterationLog(iter).objective = objectiveCurrent.total;
        iterationLog(iter).phiData = objectiveCurrent.phiData;
        iterationLog(iter).phiReg = objectiveCurrent.phiReg;
        iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
        iterationLog(iter).stepLength = best.stepLength;

        if objectiveCurrent.normalizedDataMisfit <= 1e-4
            exitReason = "targetMisfit";
            break;
        end
    end

    iterationLog = iterationLog([iterationLog.iter] > 0);
    summary = struct();
    summary.nCells = mesh.nElements;
    summary.nMeasurements = numel(dObs);
    summary.objectiveBefore = objectiveStart.total;
    summary.phiDataBefore = objectiveStart.phiData;
    summary.normalizedDataMisfitBefore = objectiveStart.normalizedDataMisfit;
    summary.rawResidualMeanAbsBefore = mean(abs(rawResidualStart));
    summary.rawResidualRmsBefore = sqrt(mean(rawResidualStart.^2));
    summary.rawResidualMaxAbsBefore = max(abs(rawResidualStart));
    summary.objectiveAfter = objectiveCurrent.total;
    summary.phiDataAfter = objectiveCurrent.phiData;
    summary.phiRegAfter = objectiveCurrent.phiReg;
    summary.normalizedDataMisfitAfter = objectiveCurrent.normalizedDataMisfit;
    summary.iterationsCompleted = numel(iterationLog);
    summary.exitReason = exitReason;
    summary.lastStepLength = lastStepLength;
    summary.lastDeltaNorm = norm(lastDelta);
    summary.firstDeltaWeightGuidance = firstDeltaGuidance;
    summary.finalAnomalyPenaltyComponents = localAnomalyPenaltyComponents(mesh, param, uAnomCurrent, uAnomRef);
    summary.trueModelMeta = [];
    summary.cellSensitivity = sqrt(sum(JuFinal.^2, 1))';
end

function [R, rhsReg] = localBuildAnomalyOnlyRegularization(mesh, param, uAnomCurrent, uAnomRef, weights)
    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    R = [ ...
        sqrt(weights.zeroMeanAnom) * param.Amean; ...
        sqrt(weights.betaAnom) * speye(numel(uAnomCurrent)); ...
        sqrt(weights.axAnom) * Dx; ...
        sqrt(weights.ayAnom) * Dy; ...
        sqrt(weights.azAnom) * Dz; ...
        sqrt(weights.cxAnom) * Dxx; ...
        sqrt(weights.cyAnom) * Dyy; ...
        sqrt(weights.czAnom) * Dzz];
    rhsReg = [ ...
        -sqrt(weights.zeroMeanAnom) * (param.Amean * uAnomCurrent); ...
        sqrt(weights.betaAnom) * (uAnomRef - uAnomCurrent); ...
        -sqrt(weights.axAnom) * (Dx * uAnomCurrent); ...
        -sqrt(weights.ayAnom) * (Dy * uAnomCurrent); ...
        -sqrt(weights.azAnom) * (Dz * uAnomCurrent); ...
        -sqrt(weights.cxAnom) * (Dxx * uAnomCurrent); ...
        -sqrt(weights.cyAnom) * (Dyy * uAnomCurrent); ...
        -sqrt(weights.czAnom) * (Dzz * uAnomCurrent)];
end

function info = localEvaluateAnomalyOnlyObjective(dObs, dPred, Wd, mesh, param, uAnom, uAnomRef, weights)
    dataResidual = Wd * (dObs - dPred);
    phiData = sum(dataResidual.^2);
    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    phiReg = weights.zeroMeanAnom * sum((param.Amean * uAnom).^2) + ...
        weights.betaAnom * sum((uAnom - uAnomRef).^2) + ...
        weights.axAnom * sum((Dx * uAnom).^2) + ...
        weights.ayAnom * sum((Dy * uAnom).^2) + ...
        weights.azAnom * sum((Dz * uAnom).^2) + ...
        weights.cxAnom * sum((Dxx * uAnom).^2) + ...
        weights.cyAnom * sum((Dyy * uAnom).^2) + ...
        weights.czAnom * sum((Dzz * uAnom).^2);
    info = struct();
    info.phiData = phiData;
    info.phiReg = phiReg;
    info.total = phiData + phiReg;
    info.normalizedDataMisfit = phiData / max(1, numel(dObs));
end

function guidance = localBuildAnomalyWeightGuidance(mesh, param, proposedAnomaly, uAnomRef, phiDataReference)
    raw = localAnomalyPenaltyComponents(mesh, param, proposedAnomaly, uAnomRef);
    names = fieldnames(raw);
    guidance = struct();
    guidance.phiDataReference = phiDataReference;
    for i = 1:numel(names)
        name = names{i};
        value = raw.(name);
        guidance.rawPenalty.(name) = value;
        if isfinite(value) && value > 0
            guidance.weightForDataScale.(name) = phiDataReference / value;
        else
            guidance.weightForDataScale.(name) = Inf;
        end
    end
end

function parts = localAnomalyPenaltyComponents(mesh, param, uAnom, uAnomRef)
    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    parts = struct();
    parts.zeroMean = sum((param.Amean * uAnom).^2);
    parts.beta = sum((uAnom - uAnomRef).^2);
    parts.ax = sum((Dx * uAnom).^2);
    parts.ay = sum((Dy * uAnom).^2);
    parts.az = sum((Dz * uAnom).^2);
    parts.cx = sum((Dxx * uAnom).^2);
    parts.cy = sum((Dyy * uAnom).^2);
    parts.cz = sum((Dzz * uAnom).^2);
end

function [mCurrent, sigmaCurrent, fieldsCurrent, iterationLog, summary] = localRunPropertyEcCalibration( ...
        mesh, sigmaTarget, mStart, mRef, weights, prior, bounds, config, aux, ecModel, maxIter, lineSearchSteps)

    mCurrent = mStart;
    [sigmaCurrent, ~, fieldsCurrent] = predict_soil_ec_physics(mCurrent, config, aux, ecModel);
    yTarget = log10(max(sigmaTarget, realmin));
    dataStd = max(0.02, 0.05 * max(abs(yTarget), 1e-3));
    Wsigma = spdiags(1 ./ dataStd, 0, numel(yTarget), numel(yTarget));

    objectiveCurrent = localEvaluatePropertyEcObjective(mCurrent, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaCurrent);
    objectiveStart = objectiveCurrent.total;
    iterationLog = repmat(struct('iter',0,'objective',0,'phiData',0,'phiReg',0,'phiTexture',0,'phiPrior',0,'normalizedDataMisfit',0,'stepNorm',0,'stepLength',0,'accepted',false), maxIter, 1);
    exitReason = "maxIter";
    lastStepLength = 0;

    for iter = 1:maxIter
        [sigmaCurrent, JsoilSigma, fieldsCurrent] = predict_soil_ec_physics(mCurrent, config, aux, ecModel);
        JlogSigma = spdiags(1 ./ (log(10) * max(sigmaCurrent, realmin)), 0, config.nCells, config.nCells) * JsoilSigma;
        [R, rhsReg] = build_property_regularization(config, mesh, mCurrent, mRef, weights, prior);
        [Rec, rhsEc, ~] = build_property_ec_prior_penalty(sigmaCurrent, JsoilSigma, prior, weights);
        R = [R; Rec];
        rhsReg = [rhsReg; rhsEc];
        rData = yTarget - log10(max(sigmaCurrent, realmin));
        [deltaM, ~] = solve_property_update(JlogSigma, rData, Wsigma, R, rhsReg);

        accepted = false;
        best = struct();
        for stepLength = lineSearchSteps
            [mTrial, fieldsTrial] = apply_property_bounds(mCurrent + stepLength * deltaM, config, bounds);
            [sigmaTrial, ~, fieldsTrial] = predict_soil_ec_physics(mTrial, config, aux, ecModel);
            objectiveTrial = localEvaluatePropertyEcObjective(mTrial, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaTrial);
            if objectiveTrial.total < objectiveCurrent.total
                accepted = true;
                best.m = mTrial;
                best.sigma = sigmaTrial;
                best.fields = fieldsTrial;
                best.objective = objectiveTrial;
                best.stepLength = stepLength;
                break;
            end
        end

        iterationLog(iter).iter = iter;
        iterationLog(iter).objective = objectiveCurrent.total;
        iterationLog(iter).phiData = objectiveCurrent.phiData;
        iterationLog(iter).phiReg = objectiveCurrent.phiReg;
        iterationLog(iter).phiTexture = objectiveCurrent.phiTexture;
        iterationLog(iter).phiPrior = objectiveCurrent.phiPrior;
        iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
        iterationLog(iter).stepNorm = norm(deltaM) / max(norm(mCurrent), eps);
        iterationLog(iter).accepted = accepted;

        if ~accepted
            exitReason = "lineSearchFailed";
            break;
        end

        mCurrent = best.m;
        sigmaCurrent = best.sigma;
        fieldsCurrent = best.fields;
        objectiveCurrent = best.objective;
        lastStepLength = best.stepLength;
        iterationLog(iter).objective = objectiveCurrent.total;
        iterationLog(iter).phiData = objectiveCurrent.phiData;
        iterationLog(iter).phiReg = objectiveCurrent.phiReg;
        iterationLog(iter).phiTexture = objectiveCurrent.phiTexture;
        iterationLog(iter).phiPrior = objectiveCurrent.phiPrior;
        iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
        iterationLog(iter).stepLength = best.stepLength;
    end

    iterationLog = iterationLog([iterationLog.iter] > 0);
    summary = struct();
    summary.exitReason = exitReason;
    summary.iterationsCompleted = numel(iterationLog);
    summary.objectiveBefore = objectiveStart;
    summary.objectiveAfter = objectiveCurrent.total;
    summary.normalizedDataMisfitAfter = objectiveCurrent.normalizedDataMisfit;
    summary.priorPenaltyAfter = objectiveCurrent.phiPrior;
    summary.priorEcPenaltyAfter = objectiveCurrent.phiPriorEC;
    summary.lastStepLength = lastStepLength;
    summary.nPriorRowsUsed = localCountPriorRows(prior);
end

function info = localEvaluatePropertyEcObjective(m, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaCurrent)
    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    fields = unpack_soil_fields(m, config);
    dataResidual = Wsigma * (log10(max(sigmaTarget, realmin)) - log10(max(sigmaCurrent, realmin)));
    phiData = sum(dataResidual.^2);

    if isfield(weights, 'depthBasis')
        basisSettings = weights.depthBasis;
    else
        basisSettings = struct('nGaussians', 3, 'includeConstant', true);
    end
    Bdepth = make_depth_basis_matrix(mesh, basisSettings.nGaussians, basisSettings.includeConstant);
    Qdepth = orth(Bdepth);
    Pperp = sparse(speye(config.nCells) - Qdepth * Qdepth');

    phiReg = 0;
    phiDepthBasis = 0;
    for i = 1:config.nProps
        name = config.names(i);
        idx = config.index.(char(name));
        u = m(idx);
        uRef = mRef(idx);
        w = weights.(char(name));
        phiReg = phiReg + w.beta * sum((u - uRef).^2);
        phiReg = phiReg + w.ax * sum((Dx * u).^2) + w.ay * sum((Dy * u).^2) + w.az * sum((Dz * u).^2);
        phiReg = phiReg + w.cx * sum((Dxx * u).^2) + w.cy * sum((Dyy * u).^2) + w.cz * sum((Dzz * u).^2);
        if isfield(w, 'depthBasisGamma') && w.depthBasisGamma > 0
            phiDepthBasis = phiDepthBasis + w.depthBasisGamma * sum((Pperp * u).^2);
        end
    end
    textureResidual = fields.Clay + fields.Sand + fields.Silt - 100;
    phiTexture = weights.textureSumGamma * sum(textureResidual.^2);

    phiPrior = 0;
    if ~isempty(prior) && isfield(prior, 'hasMeshMapping') && prior.hasMeshMapping && isfield(weights, 'prior')
        priorProps = string(fieldnames(weights.prior));
        for i = 1:numel(priorProps)
            name = priorProps(i);
            if ~ismember(name, config.names) || ~isfield(prior, char(name))
                continue;
            end
            gamma = weights.prior.(char(name));
            if ~(isnumeric(gamma) && isscalar(gamma) && gamma > 0)
                continue;
            end
            residual = fields.(char(name))(prior.cellIndex) - prior.(char(name));
            phiPrior = phiPrior + gamma * sum(residual.^2);
        end
    end

    phiPriorEC = 0;
    if ~isempty(prior) && isfield(prior, 'PointEC') && isfield(weights, 'prior') && isfield(weights.prior, 'PointEC')
        gamma = weights.prior.PointEC;
        if gamma > 0 && isfield(prior, 'cellIndex')
            ecResidual = sigmaCurrent(prior.cellIndex) - prior.PointEC;
            phiPriorEC = gamma * sum(ecResidual.^2);
        end
    end

    info = struct();
    info.total = phiData + phiReg + phiDepthBasis + phiTexture + phiPrior + phiPriorEC;
    info.phiData = phiData;
    info.phiReg = phiReg + phiDepthBasis;
    info.phiTexture = phiTexture;
    info.phiPrior = phiPrior;
    info.phiPriorEC = phiPriorEC;
    info.normalizedDataMisfit = phiData / max(1, numel(sigmaTarget));
end

function [Rbg, rhsBg] = localBuildBackgroundRegularization(param, uBgCurrent, uBgRef, weights)
    [DzBg, DzzBg] = localBackgroundDifferenceOperators(param.zCenters);
    bgProfileCurrent = param.basisZ * uBgCurrent;
    Rbg = [ ...
        sqrt(weights.betaBg) * speye(numel(uBgCurrent)); ...
        sqrt(weights.z1Bg) * (DzBg * param.basisZ); ...
        sqrt(weights.z2Bg) * (DzzBg * param.basisZ)];
    rhsBg = [ ...
        sqrt(weights.betaBg) * (uBgRef - uBgCurrent); ...
        -sqrt(weights.z1Bg) * (DzBg * bgProfileCurrent); ...
        -sqrt(weights.z2Bg) * (DzzBg * bgProfileCurrent)];
end

function phiReg = localBackgroundRegularizationValue(param, uBgCurrent, uBgRef, weights)
    [DzBg, DzzBg] = localBackgroundDifferenceOperators(param.zCenters);
    bgProfileCurrent = param.basisZ * uBgCurrent;
    phiReg = weights.betaBg * sum((uBgCurrent - uBgRef).^2) + ...
        weights.z1Bg * sum((DzBg * bgProfileCurrent).^2) + ...
        weights.z2Bg * sum((DzzBg * bgProfileCurrent).^2);
end

function [DzBg, DzzBg] = localBackgroundDifferenceOperators(zCenters)
    nZ = numel(zCenters);
    if nZ <= 1
        DzBg = sparse(0, nZ);
        DzzBg = sparse(0, nZ);
        return;
    end
    h = max(diff(zCenters(:)), eps);
    nRows1 = nZ - 1;
    rr = repelem((1:nRows1)', 2, 1);
    cc = [(1:nRows1)'; (2:nZ)'];
    vv = [-1 ./ h; 1 ./ h];
    DzBg = sparse(rr, cc, vv, nRows1, nZ);
    if nZ <= 2
        DzzBg = sparse(0, nZ);
    else
        nRows2 = nZ - 2;
        hLeft = h(1:end-1);
        hRight = h(2:end);
        denom = hLeft + hRight;
        rr = repelem((1:nRows2)', 3, 1);
        cc = [(1:nRows2)'; (2:nZ-1)'; (3:nZ)'];
        vv = [ ...
            2 ./ (hLeft .* denom); ...
            -2 .* (1 ./ hLeft + 1 ./ hRight) ./ denom; ...
            2 ./ (hRight .* denom)];
        DzzBg = sparse(rr, cc, vv, nRows2, nZ);
    end
end

function aux = localMakeAux(nCells)
    aux = struct();
    aux.pH_CaCl2 = 6.4 * ones(nCells, 1);
    aux.pH_H2O = 7.0 * ones(nCells, 1);
    aux.CaCO3 = 2.0 * ones(nCells, 1);
    aux.Coarse = 5.0 * ones(nCells, 1);
end

function surveyFile = localFindMostRecentSurveyFile(surveyDir)
    listing = dir(fullfile(surveyDir, '*.mat'));
    if isempty(listing)
        error('No electrode survey MAT files found in %s.', surveyDir);
    end
    [~, idx] = max([listing.datenum]);
    surveyFile = fullfile(listing(idx).folder, listing(idx).name);
end

function [mStart, fieldsStart] = localBuildPriorInformedStartModel(config, mesh, refFields, prior)
    fieldsStart = refFields;
    if isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'nSamples') || prior.nSamples == 0
        mStart = pack_soil_fields(fieldsStart, config);
        return;
    end
    z = mesh.elementCenters(:, 3);
    x = mesh.elementCenters(:, 1);
    y = mesh.elementCenters(:, 2);
    zBlend = 2.0;
    depthWeight = exp(-0.5 * (z ./ zBlend).^2);
    initPropNames = intersect(config.names, ["N","P","K","OC","Moisture","Clay","Sand","Silt"], 'stable');
    for i = 1:numel(initPropNames)
        name = initPropNames(i);
        if ~isfield(prior, char(name))
            continue;
        end
        base = fieldsStart.(char(name));
        sampleVals = prior.(char(name));
        Fnat = scatteredInterpolant(prior.xy(:, 1), prior.xy(:, 2), sampleVals, 'natural', 'nearest');
        surfaceVals = Fnat(x, y);
        fieldsStart.(char(name)) = base + depthWeight .* (surfaceVals - base);
    end
    if isfield(fieldsStart, 'Silt') && isfield(fieldsStart, 'Clay') && isfield(fieldsStart, 'Sand')
        fieldsStart.Silt = max(fieldsStart.Silt, 5);
        fieldsStart.Sand = max(fieldsStart.Sand, 5);
        textureSum = fieldsStart.Clay + fieldsStart.Silt + fieldsStart.Sand;
        fieldsStart.Clay = 100 * fieldsStart.Clay ./ max(textureSum, eps);
        fieldsStart.Silt = 100 * fieldsStart.Silt ./ max(textureSum, eps);
        fieldsStart.Sand = 100 * fieldsStart.Sand ./ max(textureSum, eps);
    end
    mStart = pack_soil_fields(fieldsStart, config);
end

function localRefreshAlignedPriorCsv(priorCsvPath, refFields, ecModel, truthSettings)
    if ~isfile(priorCsvPath)
        return;
    end
    T = readtable(priorCsvPath, 'TextType', 'string');
    if ~all(ismember(["x_m","y_m"], string(T.Properties.VariableNames)))
        return;
    end
    nSamples = height(T);
    points = [T.x_m, T.y_m, zeros(nSamples, 1)];
    fields = evaluate_soil_true_model_at_points(points, refFields, truthSettings);
    configPts = localMakePointPropertyConfig(nSamples);
    auxPts = localMakeAux(nSamples);
    mPts = pack_soil_fields(fields, configPts);
    [sigmaPts, ~] = predict_soil_ec_physics(mPts, configPts, auxPts, ecModel);
    T.pH_CaCl2 = auxPts.pH_CaCl2;
    T.pH_H2O = auxPts.pH_H2O;
    T.point_ec_s_m = sigmaPts;
    T.N = fields.N;
    T.P = fields.P;
    T.K = fields.K;
    T.OC = fields.OC;
    T.moisture_pct = fields.Moisture;
    T.clay_pct = fields.Clay;
    T.silt_pct = fields.Silt;
    T.sand_pct = fields.Sand;
    T.CaCO3 = auxPts.CaCO3;
    T.Coarse = auxPts.Coarse;
    writetable(T, priorCsvPath);
end

function config = localMakePointPropertyConfig(nCells)
    names = ["N","P","K","OC","Moisture","Clay","Sand","Silt"];
    nProps = numel(names);
    index = struct();
    for i = 1:nProps
        index.(char(names(i))) = ((i - 1) * nCells + 1):(i * nCells);
    end
    config = struct('names', names, 'nProps', nProps, 'nCells', nCells, ...
        'totalSize', nProps * nCells, 'index', index, ...
        'textureNames', ["Clay","Sand","Silt"], ...
        'transformedNames', ["N","P","K","OC","Moisture","Clay","Sand","Silt"]);
end

function localPrintProblemSummary(mesh, forwardMesh, survey, surveyFile, useParallel, forwardMeshType)
    sourcePairs = unique(survey.quads(:, 1:2), 'rows', 'stable');
    fprintf('\nStaged EC/property inversion summary\n');
    fprintf('  Survey file: %s\n', surveyFile);
    fprintf('  Forward mesh type: %s\n', forwardMeshType);
    fprintf('  Electrodes: %d\n', survey.nElectrodes);
    fprintf('  Measurements (ABMN): %d\n', survey.nMeasurements);
    fprintf('  Unique A-B source pairs: %d\n', size(sourcePairs, 1));
    fprintf('  Coarse inversion cells: %d\n', mesh.nElements);
    fprintf('  Fine forward cells: %d\n', forwardMesh.nElements);
    if useParallel
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('  Jacobian mode: parallel (pool not started yet)\n');
        else
            fprintf('  Jacobian mode: parallel (%d workers available)\n', pool.NumWorkers);
        end
    else
        fprintf('  Jacobian mode: serial\n');
    end
    fprintf('\n');
end

function localPrintParsedEcWeights(weights, weightsCsvPath)
    fprintf('\nParsed EC regularization weights\n');
    fprintf('  Source: %s\n', weightsCsvPath);
    fprintf('  Background: mode=%s betaBg=%.4g z1Bg_grad_per_m=%.4g z2Bg_curv_per_m2=%.4g\n', ...
        weights.backgroundMode, weights.betaBg, weights.z1Bg, weights.z2Bg);
    fprintf('  Anomaly: zeroMean=%.4g beta=%.4g ax_grad_per_m=%.4g ay_grad_per_m=%.4g az_grad_per_m=%.4g cx_curv_per_m2=%.4g cy_curv_per_m2=%.4g cz_curv_per_m2=%.4g\n\n', ...
        weights.zeroMeanAnom, weights.betaAnom, weights.axAnom, weights.ayAnom, weights.azAnom, ...
        weights.cxAnom, weights.cyAnom, weights.czAnom);
end

function localPrintAnomalyWeightGuidance(summary, weights, filePath)
    if nargin < 3
        filePath = "";
    end
    text = localFormatAnomalyWeightGuidance(summary, weights);
    fprintf('%s', text);
    if strlength(string(filePath)) > 0
        fid = fopen(filePath, 'w');
        if fid ~= -1
            cleaner = onCleanup(@() fclose(fid));
            fprintf(fid, '%s', text);
        end
    end
end

function text = localFormatAnomalyWeightGuidance(summary, weights)
    lines = strings(0, 1);
    lines(end+1) = "";
    lines(end+1) = "Stage 2 anomaly weight guidance";
    lines(end+1) = "  Stage 1 residual entering anomaly stage:";
    lines(end+1) = sprintf("    weighted phiData: %.4g", summary.phiDataBefore);
    lines(end+1) = sprintf("    weighted mean per datum: %.4g", summary.normalizedDataMisfitBefore);
    lines(end+1) = sprintf("    raw voltage mean abs: %.4g V", summary.rawResidualMeanAbsBefore);
    lines(end+1) = sprintf("    raw voltage RMS: %.4g V", summary.rawResidualRmsBefore);
    lines(end+1) = sprintf("    raw voltage max abs: %.4g V", summary.rawResidualMaxAbsBefore);

    if ~isfield(summary, 'firstDeltaWeightGuidance') || isempty(fieldnames(summary.firstDeltaWeightGuidance))
        lines(end+1) = "  No first-update guidance available.";
        lines(end+1) = "";
        text = char(strjoin(lines, newline) + newline);
        return;
    end

    names = ["zeroMean","beta","ax","ay","az","cx","cy","cz"];
    currentWeights = [weights.zeroMeanAnom, weights.betaAnom, weights.axAnom, weights.ayAnom, ...
        weights.azAnom, weights.cxAnom, weights.cyAnom, weights.czAnom];

    lines(end+1) = "  Approximate weights where each penalty matches the Stage 1 residual";
    lines(end+1) = "  for the first proposed anomaly update:";
    lines(end+1) = sprintf("    %-9s %-13s %-13s %-13s", 'term', 'current', 'matchData', 'current/match');
    for i = 1:numel(names)
        name = names(i);
        rawName = char(name);
        if isfield(summary.firstDeltaWeightGuidance.weightForDataScale, rawName)
            matchWeight = summary.firstDeltaWeightGuidance.weightForDataScale.(rawName);
        else
            matchWeight = Inf;
        end
        currentWeight = currentWeights(i);
        if isfinite(matchWeight) && matchWeight > 0
            ratio = currentWeight / matchWeight;
        elseif currentWeight == 0
            ratio = 0;
        else
            ratio = Inf;
        end
        lines(end+1) = sprintf("    %-9s %-13.4g %-13.4g %-13.4g", rawName, currentWeight, matchWeight, ratio);
    end
    lines(end+1) = "  Reading: current/match << 1 is data-dominated; ~1 is comparable; >>1 is regularization-dominated.";
    lines(end+1) = "";
    text = char(strjoin(lines, newline) + newline);
end

function meta = localStageTrueMeta(forwardMesh, truthSettings, plotSliceZ)
    c = forwardMesh.elementCenters;
    meta = struct();
    if nargin < 2 || isempty(truthSettings)
        truthSettings = struct();
    end
    if ~isfield(truthSettings, 'anomalyCenterX'), truthSettings.anomalyCenterX = mean(c(:,1)); end
    if ~isfield(truthSettings, 'anomalyCenterY'), truthSettings.anomalyCenterY = mean(c(:,2)); end
    if ~isfield(truthSettings, 'anomalyCenterZ'), truthSettings.anomalyCenterZ = 1.6; end
    meta.anomalyCenter = [truthSettings.anomalyCenterX, truthSettings.anomalyCenterY, truthSettings.anomalyCenterZ];
    meta.sliceX = meta.anomalyCenter(1);
    meta.sliceY = meta.anomalyCenter(2);
    if nargin >= 3 && ~isempty(plotSliceZ)
        meta.sliceZ = plotSliceZ;
    else
        meta.sliceZ = meta.anomalyCenter(3);
    end
end

function slices = localBuildAnalyticEcSlices(forwardMesh, refFields, truthSettings, ecModel, auxTemplate, plotSliceZ)
    meta = localStageTrueMeta(forwardMesh, truthSettings, plotSliceZ);
    xq = localMeshCenters(forwardMesh, "x");
    yq = localMeshCenters(forwardMesh, "y");
    zq = localMeshCenters(forwardMesh, "z");

    [XYx, XYy] = meshgrid(xq, yq);
    pointsXY = [XYx(:), XYy(:), meta.sliceZ * ones(numel(XYx), 1)];
    fieldsXY = evaluate_soil_true_model_at_points(pointsXY, refFields, truthSettings);
    XYv = reshape(localPredictSigmaFromFields(fieldsXY, ecModel, auxTemplate), size(XYx));

    [XZx, XZz] = meshgrid(xq, zq);
    pointsXZ = [XZx(:), meta.sliceY * ones(numel(XZx),1), XZz(:)];
    fieldsXZ = evaluate_soil_true_model_at_points(pointsXZ, refFields, truthSettings);
    XZv = reshape(localPredictSigmaFromFields(fieldsXZ, ecModel, auxTemplate), size(XZx));

    slices = struct();
    slices.xy = struct('x', XYx, 'y', XYy, 'v', XYv, 'z', meta.sliceZ, ...
        'xEdges', forwardMesh.x(:).', 'yEdges', forwardMesh.y(:).');
    slices.xz = struct('x', XZx, 'z', XZz, 'v', XZv, 'y', meta.sliceY, ...
        'xEdges', forwardMesh.x(:).', 'zEdges', forwardMesh.z(:).');
end

function centers = localMeshCenters(mesh, axisName)
    axisName = char(axisName);
    centersName = [axisName, 'Centers'];
    if isfield(mesh, centersName)
        centers = mesh.(centersName)(:).';
        return;
    end
    if isfield(mesh, axisName)
        edges = mesh.(axisName)(:).';
        centers = 0.5 * (edges(1:end-1) + edges(2:end));
        return;
    end
    error('Mesh does not define %s or %s.', axisName, centersName);
end

function n = localCountPriorRows(prior)
    if isempty(prior) || ~isfield(prior, 'nSamples')
        n = 0;
    else
        n = prior.nSamples;
    end
end

function sigma = localPredictSigmaFromFields(fields, ecModel, auxTemplate)
    nCells = numel(fields.N);
    cfg = localMakePointPropertyConfig(nCells);
    aux = struct();
    auxNames = fieldnames(auxTemplate);
    for i = 1:numel(auxNames)
        val = auxTemplate.(auxNames{i});
        if isscalar(val)
            aux.(auxNames{i}) = repmat(val, nCells, 1);
        else
            aux.(auxNames{i}) = val(1) * ones(nCells, 1);
        end
    end
    m = pack_soil_fields(fields, cfg);
    [sigma, ~] = predict_soil_ec_physics(m, cfg, aux, ecModel);
end
