clear; clc; close all;

saveArg = true;
useParallel = true;
maxIter = 3;
tolObjective = 1e-5;
tolStep = 1e-3;
targetMisfit = 1e-4;
lineSearchSteps = [8.0, 5.0, 3.0, 2.0, 1.5, 1.0, 0.75, 0.5, 0.25, 0.1];
exploratoryLineSearchSteps = [12.0, 8.0, 5.0, 3.0];
smallImprovementPatience = 3;
highMisfitEscapeThreshold = 1.0;
forwardMeshType = "fvm_octree_style"; % fvm_octree_style | fvm_rect_refined
progressMode = "text"; % text | waitbar | off
lateralBuffer = 15; % m beyond the outer electrode footprint
truthSettings = struct( ...
    'nBackground', 30, ...
    'nAnomalyAmplitude', 6, ...
    'nAnomalyFloor', 29, ...
    'clayBackground', 5, ...
    'clayLayerAmplitude', 62, ...
    'clayFloor', 5, ...
    'sandBackground', 70, ...
    'sandLayerAmplitude', 36, ...
    'sandFloor', 5);

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));
outDir = fullfile(rootDir, 'outputs', 'ert_property_inversion_demo');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

priorCsvPath = fullfile(rootDir, 'data', 'prior_samples', 'soil_property_prior_samples_aligned_dense.csv');
weightsCsvPath = fullfile(rootDir, 'matlab', 'property_regularization_weights.csv');
surveyDir = fullfile(rootDir, 'data', 'electrode_arrays');
surveyFile = localFindMostRecentSurveyFile(surveyDir);
surveyData = load(surveyFile, 'survey', 'summary');
survey = surveyData.survey;

%% Mesh and survey

surveyCenter = mean(survey.electrodes(:, 1:2), 1);
surveyRadius = max(vecnorm(survey.electrodes(:, 1:2) - surveyCenter, 2, 2));
outerRadius = surveyRadius + lateralBuffer;

sb = max(25, surveyRadius); % survey bound
or = outerRadius;
xfc = [5,10,20,sb,or];
xff = sort([2.5,xfc,mean([xfc(1:end-1);xfc(2:end)],1)]);
zfc = [0, 1, 2, 4, 6];
zff = sort([zfc,mean([zfc(1:end-1);zfc(2:end)],1)]);

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

ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

%% Property fields and priors

config = make_soil_property_config(mesh);
configForward = make_soil_property_config(forwardMesh);
[mRefBase, refFields] = make_soil_reference_model(config);
[mRefForward, refFieldsForward] = make_soil_reference_model(configForward); %#ok<ASGLU>
[mTrue, fieldsTrue] = make_soil_true_model(config, mesh, refFields, truthSettings);
[mTrueForward, fieldsTrueForward] = make_soil_true_model(configForward, forwardMesh, refFieldsForward, truthSettings);

if contains(string(priorCsvPath), "soil_property_prior_samples_aligned_dense.csv", 'IgnoreCase', true)
    localRefreshAlignedPriorCsv(priorCsvPath, refFields, ecModel, truthSettings);
end

if isfile(priorCsvPath)
    prior = load_soil_property_prior_samples(priorCsvPath, 'Mesh', mesh);
else
    warning('Prior sample file not found: %s. Continuing without prior constraints.', priorCsvPath);
    prior = [];
end

[mStart, fieldsStart] = localBuildPriorInformedStartModel(config, mesh, refFields, prior);
mRef = mStart;

aux = localMakeAux(config.nCells);
auxForward = localMakeAux(configForward.nCells);

%% True data on fine forward mesh

[sigmaTrueForward, ~, fieldsTrueForward] = predict_soil_ec_physics(mTrueForward, configForward, auxForward, ecModel);
sigmaTrue = coarseFromFine * sigmaTrueForward;

Ktrue = assemble_ert_stiffness(forwardMesh, sigmaTrueForward, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, ...
    'ProgressMode', progressMode);
resultTrue = solve_ert_forward(forwardMesh, survey, Ktrue);
dObs = resultTrue.voltage;

%% Initial model

weights = default_property_regularization_weights(weightsCsvPath);
bounds = default_property_bounds();
dataStd = max(1e-3, 0.03 * max(abs(dObs), 1e-3));
Wd = spdiags(1 ./ dataStd, 0, numel(dObs), numel(dObs));

[sigmaStart, JsoilStart, fieldsStart] = predict_soil_ec_physics(mStart, config, aux, ecModel); %#ok<ASGLU>
uSigmaStart = log10(max(sigmaStart, realmin));
sigmaForwardStart = 10 .^ (fineFromCoarse * uSigmaStart);

Kstart = assemble_ert_stiffness(forwardMesh, sigmaForwardStart, ...
    'UseParallel', useParallel, ...
    'ElementType', forwardMesh.elementType, ...
    'ShowProgress', true, ...
    'ProgressMode', progressMode);
resultStart = solve_ert_forward(forwardMesh, survey, Kstart);
dStart = resultStart.voltage;

mCurrent = mStart;
sigmaCurrent = sigmaStart;
uSigmaCurrent = uSigmaStart;
fieldsCurrent = fieldsStart;
resultCurrent = resultStart;
dPred = dStart;
objectiveCurrent = evaluate_property_inversion_objective( ...
    mCurrent, dObs, dPred, Wd, config, mesh, mRef, weights, prior, sigmaCurrent);

iterationLog = repmat(struct( ...
    'iter', 0, ...
    'objective', 0, ...
    'phiData', 0, ...
    'phiReg', 0, ...
    'phiTexture', 0, ...
    'phiPrior', 0, ...
    'normalizedDataMisfit', 0, ...
    'stepNorm', 0, ...
    'stepLength', 0, ...
    'accepted', false), maxIter, 1);

exitReason = "maxIter";
lastDeltaM = zeros(size(mCurrent));
lastStepLength = 0;
systemInfo = struct();
regInfo = struct();
J = [];
Jert = [];
Jsoil = [];
stallCount = 0;

%% Iterative inversion loop

for iter = 1:maxIter
    fprintf('Property inversion 3D iteration %d of %d\n', iter, maxIter);

    Jert = build_ert_jacobian_fd_u(forwardMesh, survey, uSigmaCurrent, fineFromCoarse, resultCurrent, ...
        'UseParallel', useParallel, ...
        'ShowProgress', true, ...
        'ProgressMode', progressMode, ...
        'RelativePerturbation', 0.02, ...
        'AbsolutePerturbation', 1e-3);

    [sigmaCurrent, JsoilSigma, fieldsCurrent] = predict_soil_ec_physics(mCurrent, config, aux, ecModel);
    sigmaScale = log(10) * max(sigmaCurrent, realmin);
    Jsoil = spdiags(1 ./ sigmaScale, 0, config.nCells, config.nCells) * JsoilSigma;
    J = Jert * Jsoil;

    [R, rhsReg, regInfo] = build_property_regularization(config, mesh, mCurrent, mRef, weights, prior);
    [Rec, rhsEc, ecPriorInfo] = build_property_ec_prior_penalty(sigmaCurrent, JsoilSigma, prior, weights);
    R = [R; Rec];
    rhsReg = [rhsReg; rhsEc];
    regInfo.nRowsPriorEC = ecPriorInfo.nRows;
    rData = dObs - dPred;
    [deltaM, systemInfo] = solve_property_update(J, rData, Wd, R, rhsReg);

    accepted = false;
    rescueAccepted = false;
    bestTrial = struct();
    bestRescueTrial = struct();
    for stepLength = lineSearchSteps
        mTrialRaw = mCurrent + stepLength * deltaM;
        [mTrial, fieldsTrial] = apply_property_bounds(mTrialRaw, config, bounds);
        [sigmaTrial, ~, fieldsTrial] = predict_soil_ec_physics(mTrial, config, aux, ecModel);
        uSigmaTrial = log10(max(sigmaTrial, realmin));
        sigmaTrialForward = 10 .^ (fineFromCoarse * uSigmaTrial);

        Ktrial = assemble_ert_stiffness(forwardMesh, sigmaTrialForward, ...
            'UseParallel', useParallel, ...
            'ElementType', forwardMesh.elementType, ...
            'ShowProgress', true, ...
            'ProgressMode', progressMode);
        resultTrial = solve_ert_forward(forwardMesh, survey, Ktrial);
        dTrial = resultTrial.voltage;
        objectiveTrial = evaluate_property_inversion_objective( ...
            mTrial, dObs, dTrial, Wd, config, mesh, mRef, weights, prior, sigmaTrial);

        if ~isfinite(objectiveTrial.total)
            continue;
        end

        if objectiveTrial.total < objectiveCurrent.total && ...
                (~accepted || objectiveTrial.total < bestTrial.objective.total)
            accepted = true;
            bestTrial.m = mTrial;
            bestTrial.sigma = sigmaTrial;
            bestTrial.uSigma = uSigmaTrial;
            bestTrial.fields = fieldsTrial;
            bestTrial.sigmaForward = sigmaTrialForward;
            bestTrial.result = resultTrial;
            bestTrial.dPred = dTrial;
            bestTrial.objective = objectiveTrial;
            bestTrial.stepLength = stepLength;
        end

        if objectiveTrial.phiData < objectiveCurrent.phiData * 0.97 && ...
                objectiveTrial.total <= objectiveCurrent.total * 2.0 && ...
                (~rescueAccepted || objectiveTrial.phiData < bestRescueTrial.objective.phiData)
            rescueAccepted = true;
            bestRescueTrial.m = mTrial;
            bestRescueTrial.sigma = sigmaTrial;
            bestRescueTrial.uSigma = uSigmaTrial;
            bestRescueTrial.fields = fieldsTrial;
            bestRescueTrial.sigmaForward = sigmaTrialForward;
            bestRescueTrial.result = resultTrial;
            bestRescueTrial.dPred = dTrial;
            bestRescueTrial.objective = objectiveTrial;
            bestRescueTrial.stepLength = stepLength;
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
    iterationLog(iter).stepLength = lastStepLength;
    iterationLog(iter).accepted = accepted;

    if ~accepted && rescueAccepted
        accepted = true;
        bestTrial = bestRescueTrial;
    end

    if ~accepted
        exitReason = "lineSearchFailed";
        lastDeltaM = deltaM;
        lastStepLength = 0;
        break;
    end

    relativeObjectiveDrop = (objectiveCurrent.total - bestTrial.objective.total) / max(objectiveCurrent.total, eps);
    relativeStepNorm = norm(bestTrial.stepLength * deltaM) / max(norm(mCurrent), eps);

    mCurrent = bestTrial.m;
    sigmaCurrent = bestTrial.sigma;
    uSigmaCurrent = bestTrial.uSigma;
    fieldsCurrent = bestTrial.fields;
    resultCurrent = bestTrial.result;
    dPred = bestTrial.dPred;
    objectiveCurrent = bestTrial.objective;
    lastDeltaM = deltaM;
    lastStepLength = bestTrial.stepLength;

    iterationLog(iter).objective = objectiveCurrent.total;
    iterationLog(iter).phiData = objectiveCurrent.phiData;
    iterationLog(iter).phiReg = objectiveCurrent.phiReg;
    iterationLog(iter).phiTexture = objectiveCurrent.phiTexture;
    iterationLog(iter).phiPrior = objectiveCurrent.phiPrior;
    iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
    iterationLog(iter).stepNorm = relativeStepNorm;
    iterationLog(iter).stepLength = bestTrial.stepLength;

    if objectiveCurrent.normalizedDataMisfit <= targetMisfit
        exitReason = "targetMisfit";
        break;
    end
    if relativeObjectiveDrop < tolObjective
        if objectiveCurrent.normalizedDataMisfit > highMisfitEscapeThreshold
            [escapeAccepted, escapeTrial] = localTryExploratoryPropertyStep( ...
                mCurrent, mRef, deltaM, exploratoryLineSearchSteps, bounds, ...
                config, aux, ecModel, fineFromCoarse, forwardMesh, survey, ...
                dObs, Wd, mesh, weights, prior, objectiveCurrent, useParallel, progressMode);
            if escapeAccepted
                mCurrent = escapeTrial.m;
                sigmaCurrent = escapeTrial.sigma;
                uSigmaCurrent = escapeTrial.uSigma;
                fieldsCurrent = escapeTrial.fields;
                resultCurrent = escapeTrial.result;
                dPred = escapeTrial.dPred;
                objectiveCurrent = escapeTrial.objective;
                lastStepLength = escapeTrial.stepLength;
                stallCount = 0;

                iterationLog(iter).objective = objectiveCurrent.total;
                iterationLog(iter).phiData = objectiveCurrent.phiData;
                iterationLog(iter).phiReg = objectiveCurrent.phiReg;
                iterationLog(iter).phiTexture = objectiveCurrent.phiTexture;
                iterationLog(iter).phiPrior = objectiveCurrent.phiPrior;
                iterationLog(iter).normalizedDataMisfit = objectiveCurrent.normalizedDataMisfit;
                iterationLog(iter).stepLength = escapeTrial.stepLength;
                continue;
            end

            stallCount = stallCount + 1;
            if stallCount < smallImprovementPatience
                continue;
            end
        end
        exitReason = "smallObjectiveImprovement";
        break;
    else
        stallCount = 0;
    end
    if relativeStepNorm < tolStep
        exitReason = "smallStep";
        break;
    end
end

iterationLog = iterationLog([iterationLog.iter] > 0);
mUpdated = mCurrent;
sigmaUpdated = sigmaCurrent;
fieldsUpdated = fieldsCurrent;
dUpdated = dPred;

%% Save and plot

summary = struct();
summary.nCells = config.nCells;
summary.nMeasurements = numel(dObs);
summary.nPropertyParams = config.totalSize;
summary.dataMisfitBefore = norm(Wd * (dObs - dStart));
summary.dataMisfitAfter = norm(Wd * (dObs - dUpdated));
summary.objectiveBefore = evaluate_property_inversion_objective( ...
    mStart, dObs, dStart, Wd, config, mesh, mRef, weights, prior, sigmaStart).total;
summary.objectiveAfter = objectiveCurrent.total;
summary.stepNorm = norm(lastStepLength * lastDeltaM) / max(norm(mUpdated), eps);
summary.lastStepLength = lastStepLength;
summary.regularizationRows = size(R, 1);
summary.textureSumResidualBefore = mean(abs(fieldsStart.Clay + fieldsStart.Sand + fieldsStart.Silt - 100));
summary.textureSumResidualAfter = mean(abs(fieldsUpdated.Clay + fieldsUpdated.Sand + fieldsUpdated.Silt - 100));
summary.iterationsCompleted = numel(iterationLog);
summary.exitReason = exitReason;
summary.normalizedDataMisfitAfter = objectiveCurrent.normalizedDataMisfit;
summary.priorPenaltyAfter = objectiveCurrent.phiPrior;
summary.priorEcPenaltyAfter = objectiveCurrent.phiPriorEC;
summary.nPriorSamples = 0;
summary.nPriorRowsUsed = 0;
summary.nPriorEcRowsUsed = 0;
summary.priorFile = '';
summary.surveyFile = surveyFile;
summary.surveySummary = surveyData.summary;
summary.forwardMeshType = forwardMeshType;
summary.progressMode = progressMode;
if ~isempty(prior)
    summary.nPriorSamples = prior.nSamples;
    summary.priorFile = prior.filePath;
end
if isfield(regInfo, 'nRowsPrior')
    summary.nPriorRowsUsed = regInfo.nRowsPrior;
end
if isfield(regInfo, 'nRowsPriorEC')
    summary.nPriorEcRowsUsed = regInfo.nRowsPriorEC;
end

fig = plot_property_inversion_demo(mesh, config, ...
    fieldsTrue, fieldsUpdated, ...
    sigmaTrue, sigmaUpdated, ...
    dObs, dUpdated, iterationLog, targetMisfit, summary, ...
    refFields, aux, ecModel, truthSettings);
figRatio = plot_property_fit_ratios(mesh, fieldsTrue, fieldsUpdated, sigmaTrue, sigmaUpdated);

if saveArg
    save(fullfile(outDir, 'ert_property_inversion_demo.mat'), ...
        'mesh', 'forwardMesh', 'survey', 'config', 'configForward', 'aux', 'auxForward', ...
        'ecModel', 'weights', 'bounds', 'prior', 'priorCsvPath', ...
        'mRefBase', 'mRef', 'mTrue', 'mTrueForward', 'mStart', 'mCurrent', 'mUpdated', ...
        'fieldsTrue', 'fieldsTrueForward', 'fieldsStart', 'fieldsUpdated', ...
        'sigmaTrue', 'sigmaTrueForward', 'sigmaStart', 'sigmaUpdated', ...
        'dObs', 'dStart', 'dUpdated', 'Jert', 'Jsoil', 'J', ...
        'R', 'rhsReg', 'regInfo', 'lastDeltaM', 'summary', 'systemInfo', ...
        'iterationLog', 'exitReason', 'fineFromCoarse', 'coarseFromFine');
    exportgraphics(fig, fullfile(outDir, 'ert_property_inversion_demo.png'), 'Resolution', 220);
    exportgraphics(figRatio, fullfile(outDir, 'ert_property_fit_ratios.png'), 'Resolution', 220);
    write_property_inversion_summary(fullfile(outDir, 'summary.txt'), summary, weights, iterationLog);
end

%% Functions

function aux = localMakeAux(nCells)
    aux = struct();
    aux.pH_CaCl2 = 6.4 * ones(nCells, 1);
    aux.pH_H2O = 7.0 * ones(nCells, 1);
    aux.CaCO3 = 2.0 * ones(nCells, 1);
    aux.Coarse = 5.0 * ones(nCells, 1);
end

function localPrintProblemSummary(mesh, forwardMesh, survey, surveyFile, useParallel, forwardMeshType)
    sourcePairs = unique(survey.quads(:, 1:2), 'rows', 'stable');
    fprintf('\n3D property inversion summary\n');
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
    fprintf('  Property parameters: %d\n', mesh.nElements * 8);
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
    depthWeight = exp(-0.5 * (z ./ zBlend) .^ 2);

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

        blend = depthWeight;
        fieldsStart.(char(name)) = base + blend .* (surfaceVals - base);
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
    if ~all(ismember(["x_m", "y_m"], string(T.Properties.VariableNames)))
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
        first = (i - 1) * nCells + 1;
        last = i * nCells;
        index.(char(names(i))) = first:last;
    end

    config = struct();
    config.names = names;
    config.nProps = nProps;
    config.nCells = nCells;
    config.totalSize = nProps * nCells;
    config.index = index;
    config.textureNames = ["Clay","Sand","Silt"];
    config.transformedNames = ["N","P","K","OC","Moisture","Clay"];
end

function [accepted, bestTrial] = localTryExploratoryPropertyStep( ...
    mCurrent, mRef, deltaM, exploratorySteps, bounds, ...
    config, aux, ecModel, fineFromCoarse, forwardMesh, survey, ...
    dObs, Wd, mesh, weights, prior, objectiveCurrent, useParallel, progressMode)

    accepted = false;
    bestTrial = struct();

    directions = {
        deltaM, ...
        mRef - mCurrent, ...
        deltaM + 0.5 * (mRef - mCurrent)
    };

    for d = 1:numel(directions)
        direction = directions{d};
        dirNorm = norm(direction);
        if ~(isfinite(dirNorm) && dirNorm > 0)
            continue;
        end

        for stepLength = exploratorySteps
            scaledDirection = direction / dirNorm;
            trialScale = stepLength * max(norm(deltaM), 1);
            mTrialRaw = mCurrent + trialScale * scaledDirection;
            [mTrial, fieldsTrial] = apply_property_bounds(mTrialRaw, config, bounds);
            [sigmaTrial, ~, fieldsTrial] = predict_soil_ec_physics(mTrial, config, aux, ecModel);
            uSigmaTrial = log10(max(sigmaTrial, realmin));
            sigmaTrialForward = 10 .^ (fineFromCoarse * uSigmaTrial);

            Ktrial = assemble_ert_stiffness(forwardMesh, sigmaTrialForward, ...
                'UseParallel', useParallel, ...
                'ElementType', forwardMesh.elementType, ...
                'ShowProgress', false, ...
                'ProgressMode', progressMode);
            resultTrial = solve_ert_forward(forwardMesh, survey, Ktrial);
            dTrial = resultTrial.voltage;
            objectiveTrial = evaluate_property_inversion_objective( ...
                mTrial, dObs, dTrial, Wd, config, mesh, mRef, weights, prior, sigmaTrial);

            if ~isfinite(objectiveTrial.total)
                continue;
            end

            if objectiveTrial.phiData < objectiveCurrent.phiData * 0.95 && ...
                    objectiveTrial.total <= objectiveCurrent.total * 4.0 && ...
                    (~accepted || objectiveTrial.phiData < bestTrial.objective.phiData)
                accepted = true;
                bestTrial.m = mTrial;
                bestTrial.sigma = sigmaTrial;
                bestTrial.uSigma = uSigmaTrial;
                bestTrial.fields = fieldsTrial;
                bestTrial.sigmaForward = sigmaTrialForward;
                bestTrial.result = resultTrial;
                bestTrial.dPred = dTrial;
                bestTrial.objective = objectiveTrial;
                bestTrial.stepLength = trialScale;
            end
        end
    end
end
