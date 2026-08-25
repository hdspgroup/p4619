clear; clc; close all;

saveArg = true;
useParallel = true;
maxIter = 4;
tolObjective = 1e-4;
tolStep = 1e-3;
targetMisfit = 1e-4;
lineSearchSteps = [1.0, 0.5, 0.25, 0.1];
parameterizationType = "background_plus_anomaly"; % background_plus_anomaly | cellwise
startModelType = "EC_uniform_high"; % EC_uniform_high | EC_uniform_low
forwardMeshType = "fvm_octree_style"; % fvm_octree_style | fvm_rect_refined
useSensitivityAlpha = false;
progressMode = "text"; % text | waitbar | off
lateralBuffer = 15; % m, beyond the outer electrode footprint
showMeshStructureFigure = false;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));

outDir = fullfile(rootDir, 'outputs', 'ert_ec_only_inversion_demo_3d');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

surveyDir = fullfile(rootDir, 'data', 'electrode_arrays');
surveyFile = localFindMostRecentSurveyFile(surveyDir);
surveyData = load(surveyFile, 'survey', 'summary');
survey = surveyData.survey;

startModelTag = regexprep(char(lower(startModelType)), '[^a-z0-9_]+', '_');

%% Mesh and survey

surveyCenter = mean(survey.electrodes(:, 1:2), 1);
surveyRadius = max(vecnorm(survey.electrodes(:, 1:2) - surveyCenter, 2, 2));
outerRadius = surveyRadius + lateralBuffer;

xfc = unique([5, 10, 20, max(25, surveyRadius), outerRadius]);
xff = unique([2.5, 5, 10, 15, max(20, surveyRadius), outerRadius]);

xCoarse = surveyCenter(1) + unique([-fliplr(xfc), 0, xfc]);
yCoarse = surveyCenter(2) + unique([-fliplr(xfc), 0, xfc]);
zCoarse = [0, 1, 2, 3, 4.5, 6, 7];

xFine = surveyCenter(1) + unique([-fliplr(xff), 0, xff]);
yFine = surveyCenter(2) + unique([-fliplr(xff), 0, xff]);
zFine = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5.5, 7];
coarseOuterRadius = outerRadius;

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
            'Radius', coarseOuterRadius, ...
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

if showMeshStructureFigure
    figMesh = plot_forward_mesh_structure_3d(forwardMesh, survey);
    if saveArg
        exportgraphics(figMesh, fullfile(outDir, ['ert_forward_mesh_structure_3d_' startModelTag '.png']), 'Resolution', 220);
    end
end

[fineFromCoarse, coarseFromFine] = make_structured_cell_mapping(mesh, forwardMesh);

localPrintProblemSummary(mesh, forwardMesh, survey, surveyFile, useParallel, forwardMeshType);

%% Synthetic true EC model on the fine forward mesh

[sigmaTrueForward, sigmaTrueForwardNodes, trueMeta] = localMakeTrueEcModel3D(forwardMesh);
sigmaTrue = coarseFromFine * sigmaTrueForward;

%% Initial EC model

sigmaStart = localBuildInitialEcModel3D(startModelType, mesh, forwardMesh, sigmaTrueForward, coarseFromFine, trueMeta);
sigmaStart = sigmaStart(:);
uStart = log10(max(sigmaStart, realmin));

if parameterizationType == "background_plus_anomaly"
    param = make_background_anomaly_parameterization(mesh, 'BackgroundMode', "gaussian3");
    [uBgRef, uAnomRef] = decompose_u_background_anomaly(mesh, uStart, param);
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

%% True data

Ktrue = assemble_ert_stiffness(forwardMesh, sigmaTrueForward, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, ...
    'ProgressMode', progressMode);
resultTrue = solve_ert_forward(forwardMesh, survey, Ktrue);
dObs = resultTrue.voltage;

%% Initial predicted data

dataStd = max(1e-3, 0.03 * max(abs(dObs), 1e-3));
Wd = spdiags(1 ./ dataStd, 0, numel(dObs), numel(dObs));

sigmaForwardCurrent = 10 .^ (fineFromCoarse * uCurrent);
Kcurrent = assemble_ert_stiffness(forwardMesh, sigmaForwardCurrent, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, ...
    'ProgressMode', progressMode);
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
    fprintf('EC-only 3D iteration %d of %d\n', iter, maxIter);

    Ju = build_ert_jacobian_fd_u(forwardMesh, survey, uCurrent, fineFromCoarse, resultCurrent, ...
        'UseParallel', useParallel, ...
        'ShowProgress', true, ...
        'ProgressMode', progressMode, ...
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
            'ShowProgress', true, ...
            'ProgressMode', progressMode);
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
summary.parameterizationType = parameterizationType;
summary.startModelType = startModelType;
summary.forwardMeshType = forwardMeshType;
summary.useSensitivityAlpha = useSensitivityAlpha;
summary.progressMode = progressMode;
summary.showMeshStructureFigure = showMeshStructureFigure;
summary.surveyFile = surveyFile;
summary.surveySummary = surveyData.summary;
summary.trueModelMeta = trueMeta;
if useSensitivityAlpha
    summary.cellSensitivity = localCellSensitivityFromJacobian(Ju, mesh);
else
    summary.cellSensitivity = [];
end

%% Plotting

fig = plot_ec_only_inversion_demo_3d( ...
    mesh, sigmaTrue, sigmaStart, sigmaUpdated, dObs, dPred, iterationLog, targetMisfit, summary, ...
    'ForwardMesh', forwardMesh, ...
    'SigmaTrueForward', sigmaTrueForward, ...
    'SigmaTrueForwardNodes', sigmaTrueForwardNodes, ...
    'TrueDisplaySlices', localBuildTrueDisplaySlices(forwardMesh, trueMeta), ...
    'Survey', survey, ...
    'TolObjective', tolObjective);

if saveArg
    exportgraphics(fig, fullfile(outDir, ['ert_ec_only_inversion_demo_3d_' startModelTag '.png']), 'Resolution', 220);
end

%% Functions

function [sigmaTrueForward, sigmaTrueForwardNodes, meta] = localMakeTrueEcModel3D(forwardMesh)
%LOCALMAKETRUEECMODEL3D Synthetic 3D model resembling the earlier 2D case.

    c = forwardMesh.elementCenters;
    sigmaTrueForward = localEvaluateTrueEcModel3D(c);
    sigmaTrueForwardNodes = localEvaluateTrueEcModel3D(forwardMesh.nodes);

    x = c(:, 1);
    y = c(:, 2);

    meta = struct();
    meta.anomalyCenter = [mean(x), mean(y), 1.6];
    meta.sliceX = meta.anomalyCenter(1);
    meta.sliceY = meta.anomalyCenter(2);
    meta.sliceZ = meta.anomalyCenter(3);
    meta.highLayerDepth = 1.9;
    meta.bedrockDepth = 3.0;
end

function sigmaVals = localEvaluateTrueEcModel3D(coords)
%LOCALEVALUATETRUEECMODEL3D Evaluate the analytic true EC model at coordinates.

    x = coords(:, 1);
    y = coords(:, 2);
    z = coords(:, 3);

    layerBoost = 0.18 * exp(-((z - 1.9).^2) / 0.7);
    bedrockDrop = 0.18 ./ (1 + exp(-(z - 3.0) / 0.18));
    background = 1.72 + layerBoost - bedrockDrop;

    anomalyCenter = [mean(x), mean(y), 1.6];
    anomalyBlob = exp(-((x - anomalyCenter(1)).^2) / 36 ...
                    - ((y - anomalyCenter(2)).^2) / 36 ...
                    - ((z - anomalyCenter(3)).^2) / 6.25);

    sigmaVals = background - 0.12 * anomalyBlob;
    sigmaVals = max(sigmaVals, 0.6);
end

function slices = localBuildTrueDisplaySlices(forwardMesh, meta)
%LOCALBUILDTRUEDISPLAYSLICES Evaluate the analytic truth on clean XY/XZ display grids.

    xMin = min(forwardMesh.nodes(:, 1));
    xMax = max(forwardMesh.nodes(:, 1));
    yMin = min(forwardMesh.nodes(:, 2));
    yMax = max(forwardMesh.nodes(:, 2));
    zMin = min(forwardMesh.nodes(:, 3));
    zMax = max(forwardMesh.nodes(:, 3));

    xq = linspace(xMin, xMax, 180);
    yq = linspace(yMin, yMax, 180);
    zq = linspace(zMin, zMax, 140);

    [XYx, XYy] = meshgrid(xq, yq);
    XYcoords = [XYx(:), XYy(:), meta.sliceZ * ones(numel(XYx), 1)];
    XYv = reshape(localEvaluateTrueEcModel3D(XYcoords), size(XYx));

    [XZx, XZz] = meshgrid(xq, zq);
    XZcoords = [XZx(:), meta.sliceY * ones(numel(XZx), 1), XZz(:)];
    XZv = reshape(localEvaluateTrueEcModel3D(XZcoords), size(XZx));

    slices = struct();
    slices.xy = struct('x', XYx, 'y', XYy, 'v', XYv, 'z', meta.sliceZ);
    slices.xz = struct('x', XZx, 'z', XZz, 'v', XZv, 'y', meta.sliceY);
end

function sigmaStart = localBuildInitialEcModel3D(startModelType, coarseMesh, forwardMesh, sigmaTrueForward, coarseFromFine, meta)
%LOCALBUILDINITIALECMODEL3D Deterministic 3D EC starting models.

    sigmaTrueCoarse = coarseFromFine * sigmaTrueForward(:);
    [highVal, lowVal] = localEstimateBackgroundLevels3D(forwardMesh, sigmaTrueForward(:), meta);
    centers = coarseMesh.elementCenters;

    switch string(startModelType)
        case "EC_uniform_high"
            sigmaStart = highVal * ones(coarseMesh.nElements, 1);

        case "EC_uniform_low"
            sigmaStart = lowVal * ones(coarseMesh.nElements, 1);

        otherwise
            error('Unsupported startModelType: %s', startModelType);
    end

    sigmaStart = max(min(sigmaStart, max(sigmaTrueCoarse) * 5), min(sigmaTrueCoarse) * 0.1);
    sigmaStart = sigmaStart(:);
end

function [highVal, lowVal] = localEstimateBackgroundLevels3D(forwardMesh, sigmaTrueForward, meta)
%LOCALESTIMATEBACKGROUNDLEVELS3D Estimate conductive-layer and bedrock background levels.

    c = forwardMesh.elementCenters;
    x = c(:, 1);
    y = c(:, 2);
    z = c(:, 3);

    radialDist = hypot(x - meta.anomalyCenter(1), y - meta.anomalyCenter(2));
    backgroundMask = radialDist >= 35;
    if ~any(backgroundMask)
        backgroundMask = true(size(x));
    end

    highMask = backgroundMask & z >= 1.0 & z <= 3.0;
    if any(highMask)
        highVal = mean(sigmaTrueForward(highMask));
    else
        highVal = mean(sigmaTrueForward(backgroundMask));
    end

    deepMask = backgroundMask & z >= 3.5;
    if any(deepMask)
        lowVal = mean(sigmaTrueForward(deepMask));
    else
        lowVal = min(sigmaTrueForward(backgroundMask));
    end
end

function localPrintProblemSummary(mesh, forwardMesh, survey, surveyFile, useParallel, forwardMeshType)
%LOCALPRINTPROBLEMSUMMARY Print the main dimensions driving 3D runtime.

    sourcePairs = unique(survey.quads(:, 1:2), 'rows', 'stable');
    fprintf('\n3D ERT problem summary\n');
    fprintf('  Survey file: %s\n', surveyFile);
    fprintf('  Forward mesh type: %s\n', forwardMeshType);
    fprintf('  Electrodes: %d\n', survey.nElectrodes);
    fprintf('  Measurements (ABMN): %d\n', survey.nMeasurements);
    fprintf('  Unique A-B source pairs: %d\n', size(sourcePairs, 1));
    fprintf('  Coarse inversion cells: %d\n', mesh.nElements);
    fprintf('  Fine forward cells: %d\n', forwardMesh.nElements);
    fprintf('  Coarse mesh size [nx ny nz]: [%d %d %d]\n', mesh.gridSize(1) - 1, mesh.gridSize(2) - 1, mesh.gridSize(3) - 1);
    if isfield(forwardMesh, 'gridSize')
        fprintf('  Fine mesh size   [nx ny nz]: [%d %d %d]\n', ...
            forwardMesh.gridSize(1) - 1, forwardMesh.gridSize(2) - 1, forwardMesh.gridSize(3) - 1);
    else
        fprintf('  Fine mesh type: unstructured %s (%d nodes, %d elements)\n', ...
            forwardMesh.elementType, forwardMesh.nNodes, forwardMesh.nElements);
    end
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

function surveyFile = localFindMostRecentSurveyFile(surveyDir)
%LOCALFINDMOSTRECENTSURVEYFILE Find the newest saved electrode survey MAT file.

    listing = dir(fullfile(surveyDir, '*.mat'));
    if isempty(listing)
        error('No electrode survey MAT files found in %s.', surveyDir);
    end

    [~, idx] = max([listing.datenum]);
    surveyFile = fullfile(listing(idx).folder, listing(idx).name);
end

function sens = localCellSensitivityFromJacobian(Ju, mesh)
%LOCALCELLSENSITIVITYFROMJACOBIAN Rough per-cell sensitivity for plotting alpha.

    if isempty(Ju)
        sens = ones(mesh.nElements, 1);
        return;
    end

    sens = sqrt(sum(Ju .^ 2, 1))';
    sens(~isfinite(sens)) = 0;
end
