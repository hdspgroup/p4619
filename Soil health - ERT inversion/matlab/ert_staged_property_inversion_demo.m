clear; clc; close all;

saveArg = true;
maxIterPropertyBackground = 5;
maxIterPropertyAnomaly = 10;
runPropertyBackgroundInversion = false; % false keeps Stage 1 as the prior-based background model
searchspace = 10.^(-3:0.25:3);
lineSearchPropertyBackground = searchspace;
lineSearchPropertyAnomaly = searchspace;
% lineSearchPropertyBackground = [10, 4, 2, 1, 1/2, 1/4, 1/10];
% lineSearchPropertyAnomaly = [10, 4, 2, 1, 1/2, 1/4, 1/10];
ecTargetMode = "actual"; % inversion | actual
activeZoneMode = "electrodeDistance"; % all | electrodeDistance
activeZoneDistanceMeters = []; % [] uses max A-B source spacing
propertyCropMode = "activeBoundingBox"; % activeBoundingBox | none
propertyCropPaddingCells = [0 0 0]; % [x y z] extra cells around active block
propertyPlotMaxScaleFactor = 10; % Inf uses full actual+current color range; finite values cap range vs actual span
propertyResidualOverviewScaleMode = "shared"; % shared | perProperty
propertyResidualColorLimit = 10^-1*[-1,1]; % [] auto signed-log; scalar 1 or vector [-1 1] sets direct residual CLim
hardPriorSamplesInAnomalyStage = true; % true pins sampled cells exactly to prior CSV values during stage 2
propertyAnomalyUpdateNames = ["N"]; % properties allowed to change in Stage 2; use "all" to allow every property
propertyAnomalyReferenceMode = "priorStart"; % priorStart | background

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(scriptDir, 'scr')));

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
ecBundleFile = fullfile(outDir, 'staged_ec_inversion_output.mat');
if ~isfile(ecBundleFile)
    ecBundleFile = fullfile(outDir, 'ert_staged_ec_property_inversion_demo.mat');
end
if ~isfile(ecBundleFile)
    error('Run ert_staged_ec_inversion_demo first. Missing EC bundle in %s.', outDir);
end

S = load(ecBundleFile);

priorCsvPath = fullfile(rootDir, 'data', 'prior_samples', 'soil_property_prior_samples_aligned_dense.csv');
weightsCsvPath = fullfile(rootDir, 'matlab', 'property_regularization_weights.csv');
ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

[meshFull, sigmaTargetFull, targetFieldsTrue, ecTargetMeshSource] = localSelectEcTarget(S, ecTargetMode);
activeMaskFull = localBuildActiveSoilMask(meshFull, S.survey, activeZoneMode, activeZoneDistanceMeters);
[mesh, sigmaTarget, activeMask, cropInfo] = localCropPropertyMesh( ...
    meshFull, sigmaTargetFull, activeMaskFull, propertyCropMode, propertyCropPaddingCells);
config = make_soil_property_config(mesh);
[mRefHomogeneous, refFields] = make_soil_reference_model(config);
aux = localMakeAux(config.nCells);
localVerifyPropertyGridConsistency(mesh, config, sigmaTarget, activeMask, cropInfo);
localPrintActiveSoilMask(activeMask, mesh, S.survey, activeZoneMode, activeZoneDistanceMeters, ecTargetMode, cropInfo, ecTargetMeshSource);

if contains(string(priorCsvPath), "soil_property_prior_samples_aligned_dense.csv", 'IgnoreCase', true)
    localRefreshAlignedPriorCsv(priorCsvPath, refFields, ecModel, S.truthSettings);
end

prior = load_soil_property_prior_samples(priorCsvPath, 'Mesh', mesh);
weights = default_property_regularization_weights(weightsCsvPath);
localPrintParsedPropertyWeights(weights, weightsCsvPath);
bounds = default_property_bounds();

[mStart, fieldsStart] = localBuildPriorInformedStartModel(config, mesh, refFields, prior);

fprintf('\nStaged soil-property inversion\n');
fprintf('  EC target bundle: %s\n', ecBundleFile);
fprintf('  Prior samples: %s\n', priorCsvPath);
fprintf('  Weights: %s\n', weightsCsvPath);
fprintf('  Stage 1: property background\n');

Pbg = localBuildPropertyProjection(mesh, config, weights, "background");
fprintf('  Property background basis: %d Gaussian(s), includeConstant=%d, coefficients/property=%d\n', ...
    weights.depthBasis.nGaussians, double(weights.depthBasis.includeConstant), ...
    size(Pbg.matrix, 2) / config.nProps);
if localIsConstantBackgroundProjection(Pbg, config)
    [mStartBackground, fieldsStartBackground] = localBuildHomogeneousPriorMeanStartModel(config, refFields, prior);
    fprintf('  Constant property background start: using prior mean values for each property.\n');
else
    mStartBackground = localProjectModelToProjection(mStart, Pbg);
    fieldsStartBackground = unpack_soil_fields(mStartBackground, config);
end
localPrintPropertyLayeringReport("projected property start", mesh, fieldsStartBackground, config);
weightsBg = localPropertyStageWeights(weights, "background");
if runPropertyBackgroundInversion
    [mPropertyBg, sigmaPropertyBg, fieldsPropertyBg, iterationLogPropertyBg, propertyBgSummary] = ...
        localRunPropertyEcCalibrationProjected(mesh, sigmaTarget, mStartBackground, mStartBackground, weightsBg, prior, bounds, ...
        config, aux, ecModel, Pbg, activeMask, maxIterPropertyBackground, lineSearchPropertyBackground, "property background", false);
else
    fprintf('  Property background inversion skipped: using prior-based background model.\n');
    mPropertyBg = mStartBackground;
    fieldsPropertyBg = fieldsStartBackground;
    [sigmaPropertyBg, ~, fieldsPropertyBg] = predict_soil_ec_physics(mPropertyBg, config, aux, ecModel);
    iterationLogPropertyBg = repmat(localEmptyIterationRecord(), 0, 1);
    propertyBgSummary = localSkippedPropertySummary("skippedPriorBackground", sigmaPropertyBg, sigmaTarget, activeMask);
end
propertyBgFieldHealth = localPrintFieldHealthReport("property background", fieldsPropertyBg, sigmaPropertyBg, activeMask, config, bounds);
localPrintPropertyLayeringReport("property background", mesh, fieldsPropertyBg, config);

fprintf('  Stage 2: anomaly offsets around frozen property background\n');
Panom = localBuildPropertyProjection(mesh, config, weights, "anomaly", propertyAnomalyUpdateNames);
fprintf('  Property anomaly projection: %s\n', Panom.mode);
fprintf('  Stage 2 updated properties: %s\n', strjoin(string(Panom.updateNames), ', '));
if isfield(Panom, 'matrix')
    fprintf('  Stage 2 solved variables: %d / %d\n', size(Panom.matrix, 2), config.totalSize);
end
switch lower(string(propertyAnomalyReferenceMode))
    case "priorstart"
        mPropertyAnomalyRef = mStart;
    case "background"
        mPropertyAnomalyRef = mPropertyBg;
    otherwise
        error('Unsupported propertyAnomalyReferenceMode: %s', propertyAnomalyReferenceMode);
end
fprintf('  Stage 2 regularization reference: %s\n', propertyAnomalyReferenceMode);
fprintf('  Hard prior sample pinning in anomaly stage: %d\n', double(hardPriorSamplesInAnomalyStage));
weightsAnom = localPropertyStageWeights(weights, "anomaly");
[mProp, sigmaProp, fieldsProp, iterationLogProp, propSummary] = ...
    localRunPropertyEcCalibrationProjected(mesh, sigmaTarget, mPropertyBg, mPropertyAnomalyRef, weightsAnom, prior, bounds, ...
    config, aux, ecModel, Panom, activeMask, maxIterPropertyAnomaly, lineSearchPropertyAnomaly, "property anomaly", hardPriorSamplesInAnomalyStage);
propertyFieldHealth = localPrintFieldHealthReport("property final", fieldsProp, sigmaProp, activeMask, config, bounds);
localPrintPropertyLayeringReport("property final", mesh, fieldsProp, config);
iterationLogPropertyCombined = localCombineIterationLogs(iterationLogPropertyBg, iterationLogProp);
propertyResidualReport = localPrintPropertyResidualContributions( ...
    mesh, config, fieldsProp, refFields, S.truthSettings, activeMask, weightsAnom);
priorFitReport = localPrintPriorFitReport(config, fieldsProp, prior, weights);

figProp = plot_staged_property_calibration_demo( ...
    mesh, config, targetFieldsTrue, fieldsProp, sigmaTarget, sigmaProp, ...
    iterationLogPropertyCombined, propSummary, refFields, aux, ecModel, S.truthSettings, ...
    activeMask, propertyPlotMaxScaleFactor, propertyResidualOverviewScaleMode, propertyResidualColorLimit);

if saveArg
    exportgraphics(figProp.npkOc, fullfile(outDir, 'ert_staged_property_npk_oc_demo.png'), 'Resolution', 220);
    exportgraphics(figProp.textureMoisture, fullfile(outDir, 'ert_staged_property_texture_moisture_demo.png'), 'Resolution', 220);
    exportgraphics(figProp.residualOverview, fullfile(outDir, 'ert_staged_property_residual_overview_demo.png'), 'Resolution', 220);
    exportgraphics(figProp.diagnostics, fullfile(outDir, 'ert_staged_property_diagnostics_demo.png'), 'Resolution', 220);
    exportgraphics(figProp.npkOc, fullfile(outDir, 'ert_staged_property_calibration_demo.png'), 'Resolution', 220);
    save(fullfile(outDir, 'staged_property_inversion_output.mat'), ...
        'mesh', 'config', 'mStart', 'fieldsStart', 'mStartBackground', 'fieldsStartBackground', ...
        'mRefHomogeneous', 'mPropertyBg', 'fieldsPropertyBg', ...
        'sigmaPropertyBg', 'mProp', 'fieldsProp', 'sigmaProp', ...
        'sigmaTarget', 'ecTargetMode', 'ecTargetMeshSource', 'activeZoneMode', 'activeZoneDistanceMeters', ...
        'runPropertyBackgroundInversion', ...
        'activeMask', 'propertyCropMode', 'propertyCropPaddingCells', ...
        'propertyPlotMaxScaleFactor', 'propertyResidualOverviewScaleMode', 'propertyResidualColorLimit', ...
        'hardPriorSamplesInAnomalyStage', 'propertyAnomalyUpdateNames', 'propertyAnomalyReferenceMode', 'cropInfo', ...
        'iterationLogPropertyBg', 'iterationLogProp', 'iterationLogPropertyCombined', ...
        'propertyBgSummary', 'propSummary', 'propertyResidualReport', ...
        'propertyBgFieldHealth', 'propertyFieldHealth', ...
        'propertyResidualReport', 'priorFitReport', ...
        'priorCsvPath', 'weightsCsvPath');
end

fprintf('Done. Property outputs written to %s\n', outDir);

function P = localBuildPropertyProjection(mesh, config, weights, mode, updateNames)
    if nargin < 5 || isempty(updateNames)
        updateNames = "all";
    end
    basisSettings = weights.depthBasis;
    B = make_depth_basis_matrix(mesh, basisSettings.nGaussians, basisSettings.includeConstant);
    Q = orth(B);
    switch lower(string(mode))
        case "background"
            Pcell = sparse(Q);
            blocks = cell(config.nProps, 1);
            for i = 1:config.nProps
                blocks{i} = Pcell;
            end
            P = struct('mode', "matrix", 'matrix', blkdiag(blocks{:}), ...
                'nCells', config.nCells, 'nProps', config.nProps);
        case "anomaly"
            % Solve directly in the allowed anomaly-update subspace. This is
            % important: the data, prior, and regularization rows must all see
            % the same variables that will actually be applied to the model.
            % Post-processing a full-space update after the solve makes the
            % regularization weights misleading.
            [updateMask, resolvedNames] = localPropertyUpdateMask(config, updateNames);
            cols = find(updateMask);
            Psel = sparse(cols, 1:numel(cols), 1, config.totalSize, numel(cols));
            P = struct('mode', "matrix", 'matrix', Psel, ...
                'nCells', config.nCells, 'nProps', config.nProps, ...
                'gridSize', mesh.gridSize - 1, ...
                'updateMask', updateMask, 'updateNames', resolvedNames);
        otherwise
            error('Unsupported property projection mode: %s', mode);
    end
end

function mProjected = localProjectModelToProjection(m, projection)
    mode = string(projection.mode);
    switch mode
        case "matrix"
            P = projection.matrix;
            mProjected = P * (P' * m);
        otherwise
            mProjected = m;
    end
end

function tf = localIsConstantBackgroundProjection(projection, config)
    tf = string(projection.mode) == "matrix" && ...
        isfield(projection, 'matrix') && ...
        size(projection.matrix, 2) == config.nProps;
end

function [mStart, fieldsStart] = localBuildHomogeneousPriorMeanStartModel(config, refFields, prior)
    fieldsStart = refFields;
    for name = config.names
        key = char(name);
        if isstruct(prior) && isfield(prior, key) && ~isempty(prior.(key))
            vals = prior.(key);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                fieldsStart.(key) = mean(vals) * ones(config.nCells, 1);
            end
        end
    end

    if isfield(fieldsStart, 'Clay') && isfield(fieldsStart, 'Silt') && isfield(fieldsStart, 'Sand')
        fieldsStart.Clay = max(fieldsStart.Clay, 1);
        fieldsStart.Silt = max(fieldsStart.Silt, 1);
        fieldsStart.Sand = max(fieldsStart.Sand, 1);
        textureSum = fieldsStart.Clay + fieldsStart.Silt + fieldsStart.Sand;
        fieldsStart.Clay = 100 * fieldsStart.Clay ./ max(textureSum, eps);
        fieldsStart.Silt = 100 * fieldsStart.Silt ./ max(textureSum, eps);
        fieldsStart.Sand = 100 * fieldsStart.Sand ./ max(textureSum, eps);
    end

    mStart = pack_soil_fields(fieldsStart, config);
end

function stageWeights = localPropertyStageWeights(weights, mode)
    stageWeights = weights;
    propNames = setdiff(string(fieldnames(weights)), ["depthBasis","textureSumGamma","prior"], 'stable');
    for i = 1:numel(propNames)
        name = char(propNames(i));
        w = weights.(name);
        switch lower(string(mode))
            case "background"
                stageWeights.(name).beta = localGetField(w, 'betaBg', w.beta);
                stageWeights.(name).ax = 0;
                stageWeights.(name).ay = 0;
                stageWeights.(name).az = localGetField(w, 'z1Bg', w.az);
                stageWeights.(name).cx = 0;
                stageWeights.(name).cy = 0;
                stageWeights.(name).cz = localGetField(w, 'z2Bg', w.cz);
                stageWeights.(name).depthBasisGamma = 0;
            case "anomaly"
                stageWeights.(name).depthBasisGamma = 0;
        end
    end
end

function value = localGetField(s, name, defaultValue)
    if isfield(s, name)
        value = s.(name);
    else
        value = defaultValue;
    end
end

function combined = localCombineIterationLogs(bgLog, anomalyLog)
    if isempty(bgLog)
        combined = anomalyLog(:);
        return;
    end
    combined = bgLog(:);
    if isempty(anomalyLog)
        return;
    end
    offset = numel(combined);
    anomalyLog = anomalyLog(:);
    for i = 1:numel(anomalyLog)
        anomalyLog(i).iter = offset + i;
    end
    combined = [combined; anomalyLog];
end

function rec = localEmptyIterationRecord()
    rec = struct('iter',0,'objective',0,'phiData',0,'phiReg',0,'phiTexture',0, ...
        'phiPrior',0,'normalizedDataMisfit',0,'stepNorm',0,'stepLength',0,'accepted',false);
end

function summary = localSkippedPropertySummary(exitReason, sigmaCurrent, sigmaTarget, activeMask)
    yCurrent = log10(max(sigmaCurrent(activeMask), realmin));
    yTarget = log10(max(sigmaTarget(activeMask), realmin));
    residual = yTarget - yCurrent;
    summary = struct();
    summary.exitReason = string(exitReason);
    summary.iterationsCompleted = 0;
    summary.objectiveBefore = sum(residual .^ 2);
    summary.objectiveAfter = summary.objectiveBefore;
    summary.normalizedDataMisfitAfter = summary.objectiveAfter / max(1, nnz(activeMask));
    summary.priorPenaltyAfter = 0;
    summary.priorEcPenaltyAfter = 0;
    summary.lastStepLength = 0;
    summary.priorRows = 0;
    summary.nPriorRowsUsed = 0;
end

function [meshTarget, sigmaTarget, fieldsTarget, meshSource] = localSelectEcTarget(S, ecTargetMode)
    switch lower(string(ecTargetMode))
        case "inversion"
            meshTarget = S.mesh;
            sigmaTarget = S.sigmaStage2;
            meshSource = "coarse inversion mesh";
            if isfield(S, 'fieldsTrue')
                fieldsTarget = S.fieldsTrue;
            else
                fieldsTarget = struct();
            end
        case "actual"
            if isfield(S, 'forwardMesh') && isfield(S, 'sigmaTrueForward')
                meshTarget = S.forwardMesh;
                sigmaTarget = S.sigmaTrueForward;
                meshSource = "fine forward mesh";
                if isfield(S, 'fieldsTrueForward')
                    fieldsTarget = S.fieldsTrueForward;
                elseif isfield(S, 'fieldsTrue')
                    fieldsTarget = S.fieldsTrue;
                else
                    fieldsTarget = struct();
                end
            elseif isfield(S, 'sigmaTrue')
                meshTarget = S.mesh;
                sigmaTarget = S.sigmaTrue;
                meshSource = "coarse inversion mesh (fallback; sigmaTrueForward not found)";
                if isfield(S, 'fieldsTrue')
                    fieldsTarget = S.fieldsTrue;
                else
                    fieldsTarget = struct();
                end
            else
                error('EC target mode "actual" requires sigmaTrueForward or sigmaTrue in the EC output bundle.');
            end
        otherwise
            error('Unsupported ecTargetMode: %s', ecTargetMode);
    end
    sigmaTarget = sigmaTarget(:);
    if numel(sigmaTarget) ~= meshTarget.nElements
        error('Selected EC target has %d values but selected target mesh has %d elements.', ...
            numel(sigmaTarget), meshTarget.nElements);
    end
end

function activeMask = localBuildActiveSoilMask(mesh, survey, activeZoneMode, activeZoneDistanceMeters)
    switch lower(string(activeZoneMode))
        case "all"
            activeMask = true(mesh.nElements, 1);
            return;
        case "electrodedistance"
            if isempty(activeZoneDistanceMeters)
                activeZoneDistanceMeters = localMaxSourceSpacing(survey);
            end
            electrodePts = [survey.electrodes(:, 1:2), zeros(survey.nElectrodes, 1)];
            c = mesh.elementCenters;
            minDist = inf(mesh.nElements, 1);
            for i = 1:size(electrodePts, 1)
                d = sqrt(sum((c - electrodePts(i, :)).^2, 2));
                minDist = min(minDist, d);
            end
            activeMask = minDist <= activeZoneDistanceMeters;
        otherwise
            error('Unsupported activeZoneMode: %s', activeZoneMode);
    end
end

function [meshCrop, sigmaCrop, activeCrop, cropInfo] = localCropPropertyMesh(mesh, sigmaTarget, activeMask, cropMode, paddingCells)
    cropInfo = struct();
    cropInfo.mode = char(string(cropMode));
    cropInfo.originalGridSize = mesh.gridSize;
    cropInfo.originalNElements = mesh.nElements;

    if nargin < 5 || isempty(paddingCells)
        paddingCells = [0 0 0];
    end
    paddingCells = round(double(paddingCells(:)'));
    if numel(paddingCells) == 1
        paddingCells = repmat(paddingCells, 1, 3);
    end
    paddingCells = max(paddingCells(1:3), 0);

    switch lower(string(cropMode))
        case "none"
            meshCrop = mesh;
            sigmaCrop = sigmaTarget(:);
            activeCrop = activeMask(:);
            keep = true(mesh.nElements, 1);
            cropInfo.keepCellIndex = find(keep);
            cropInfo.croppedNElements = mesh.nElements;
            cropInfo.croppedGridSize = mesh.gridSize;
            cropInfo.ixRange = [1, mesh.gridSize(1) - 1];
            cropInfo.iyRange = [1, mesh.gridSize(2) - 1];
            cropInfo.izRange = [1, mesh.gridSize(3) - 1];
            return;
        case "activeboundingbox"
            n = mesh.gridSize - 1;
            mask3 = reshape(activeMask(:), n);
            [ix, iy, iz] = ind2sub(n, find(mask3));
            if isempty(ix)
                error('Cannot crop property mesh because activeMask has no active cells.');
            end
            ixRange = [max(1, min(ix) - paddingCells(1)), min(n(1), max(ix) + paddingCells(1))];
            iyRange = [max(1, min(iy) - paddingCells(2)), min(n(2), max(iy) + paddingCells(2))];
            izRange = [max(1, min(iz) - paddingCells(3)), min(n(3), max(iz) + paddingCells(3))];

            keep3 = false(n);
            keep3(ixRange(1):ixRange(2), iyRange(1):iyRange(2), izRange(1):izRange(2)) = true;
            keep = keep3(:);
            meshCrop = localStructuredSubmesh(mesh, ixRange, iyRange, izRange);
            sigmaCrop = sigmaTarget(keep);
            activeCrop = activeMask(keep);

            cropInfo.keepCellIndex = find(keep);
            cropInfo.croppedNElements = meshCrop.nElements;
            cropInfo.croppedGridSize = meshCrop.gridSize;
            cropInfo.ixRange = ixRange;
            cropInfo.iyRange = iyRange;
            cropInfo.izRange = izRange;
        otherwise
            error('Unsupported propertyCropMode: %s', cropMode);
    end
end

function meshOut = localStructuredSubmesh(mesh, ixRange, iyRange, izRange)
    x = mesh.x(ixRange(1):ixRange(2) + 1);
    y = mesh.y(iyRange(1):iyRange(2) + 1);
    z = mesh.z(izRange(1):izRange(2) + 1);
    xc = 0.5 * (x(1:end-1) + x(2:end));
    yc = 0.5 * (y(1:end-1) + y(2:end));
    zc = 0.5 * (z(1:end-1) + z(2:end));
    [Xc, Yc, Zc] = ndgrid(xc, yc, zc);

    meshOut = struct();
    meshOut.elementType = localCopyField(mesh, 'elementType', "fvm_rect");
    meshOut.discretization = localCopyField(mesh, 'discretization', "fvm");
    meshOut.meshStyle = localCopyField(mesh, 'meshStyle', "rect_refined_crop");
    meshOut.nodes = [];
    meshOut.elements = [];
    meshOut.elementCenters = [Xc(:), Yc(:), Zc(:)];
    meshOut.x = x;
    meshOut.y = y;
    meshOut.z = z;
    meshOut.xCenters = xc(:);
    meshOut.yCenters = yc(:);
    meshOut.zCenters = zc(:);
    meshOut.dx = diff(x(:)');
    meshOut.dy = diff(y(:)');
    meshOut.dz = diff(z(:)');
    meshOut.gridSize = [numel(x), numel(y), numel(z)];
    meshOut.nCellsXYZ = meshOut.gridSize - 1;
    meshOut.nNodes = 0;
    meshOut.nElements = prod(meshOut.nCellsXYZ);
    meshOut.nCells = meshOut.nElements;
    [DX, DY, DZ] = ndgrid(meshOut.dx, meshOut.dy, meshOut.dz);
    meshOut.cellVolumes = (DX .* DY .* DZ);
    meshOut.cellVolumes = meshOut.cellVolumes(:);
end

function value = localCopyField(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function dMax = localMaxSourceSpacing(survey)
    sourcePairs = unique(survey.quads(:, 1:2), 'rows');
    a = survey.electrodes(sourcePairs(:, 1), 1:2);
    b = survey.electrodes(sourcePairs(:, 2), 1:2);
    dMax = max(vecnorm(a - b, 2, 2));
end

function localPrintActiveSoilMask(activeMask, mesh, survey, activeZoneMode, activeZoneDistanceMeters, ecTargetMode, cropInfo, ecTargetMeshSource)
    if isempty(activeZoneDistanceMeters) && strcmpi(string(activeZoneMode), "electrodeDistance")
        activeZoneDistanceMeters = localMaxSourceSpacing(survey);
    end
    fprintf('\nSoil-property EC target and active zone\n');
    fprintf('  EC target mode: %s\n', ecTargetMode);
    if nargin >= 8 && strlength(string(ecTargetMeshSource)) > 0
        fprintf('  EC target mesh: %s\n', ecTargetMeshSource);
    end
    fprintf('  Active zone mode: %s\n', activeZoneMode);
    if strcmpi(string(activeZoneMode), "electrodeDistance")
        fprintf('  Active distance: %.4g m\n', activeZoneDistanceMeters);
    end
    if nargin >= 7 && ~isempty(cropInfo)
        fprintf('  Property crop mode: %s\n', cropInfo.mode);
        fprintf('  Property cells: %d / %d retained (%.1f%% of full EC mesh)\n', ...
            cropInfo.croppedNElements, cropInfo.originalNElements, ...
            100 * cropInfo.croppedNElements / cropInfo.originalNElements);
        fprintf('  Crop cell ranges [x y z]: [%d:%d] [%d:%d] [%d:%d]\n', ...
            cropInfo.ixRange(1), cropInfo.ixRange(2), cropInfo.iyRange(1), cropInfo.iyRange(2), ...
            cropInfo.izRange(1), cropInfo.izRange(2));
    end
    fprintf('  Active cells in property mesh: %d / %d (%.1f%%)\n\n', ...
        nnz(activeMask), mesh.nElements, 100 * nnz(activeMask) / mesh.nElements);
end

function localVerifyPropertyGridConsistency(mesh, config, sigmaTarget, activeMask, cropInfo)
    cellGrid = mesh.gridSize - 1;
    nGrid = prod(cellGrid);
    if nGrid ~= mesh.nElements
        error('Property mesh grid/product mismatch: prod(gridSize-1)=%d but nElements=%d.', nGrid, mesh.nElements);
    end
    if config.nCells ~= mesh.nElements
        error('Soil config/grid mismatch: config.nCells=%d but mesh.nElements=%d.', config.nCells, mesh.nElements);
    end
    if numel(sigmaTarget) ~= config.nCells
        error('EC target/property grid mismatch: numel(sigmaTarget)=%d but config.nCells=%d.', numel(sigmaTarget), config.nCells);
    end
    if numel(activeMask) ~= config.nCells
        error('Active mask/property grid mismatch: numel(activeMask)=%d but config.nCells=%d.', numel(activeMask), config.nCells);
    end

    fprintf('\nSoil-property grid consistency check\n');
    fprintf('  Property cell grid [nx ny nz]: [%d %d %d]\n', cellGrid(1), cellGrid(2), cellGrid(3));
    fprintf('  EC target values: %d\n', numel(sigmaTarget));
    fprintf('  Soil cells per property: %d\n', config.nCells);
    fprintf('  Stacked soil unknowns: %d properties x %d cells = %d\n', ...
        config.nProps, config.nCells, config.totalSize);
    if nargin >= 5 && ~isempty(cropInfo)
        fullGrid = cropInfo.originalGridSize - 1;
        fprintf('  Full EC cell grid before crop [nx ny nz]: [%d %d %d]\n', ...
            fullGrid(1), fullGrid(2), fullGrid(3));
    end
end

function localPrintParsedPropertyWeights(weights, weightsCsvPath)
    fprintf('\nParsed soil-property regularization weights\n');
    fprintf('  Source: %s\n', weightsCsvPath);
    fprintf('  Globals: textureSumGamma=%.4g priorPointEC=%.4g depthBasisNGaussians=%g includeConstant=%g\n', ...
        weights.textureSumGamma, weights.prior.PointEC, weights.depthBasis.nGaussians, ...
        double(weights.depthBasis.includeConstant));
    fprintf('  %-10s %9s %14s %16s %9s %14s %14s %14s %9s\n', ...
        'property', 'betaBg', 'z1Bg/m', 'z2Bg/m2', 'beta', 'ax/m', 'ay/m', 'az/m', 'prior');
    propNames = ["N","P","K","OC","Moisture","Clay","Sand","Silt"];
    for i = 1:numel(propNames)
        name = char(propNames(i));
        w = weights.(name);
        priorGamma = NaN;
        if isfield(weights.prior, name)
            priorGamma = weights.prior.(name);
        end
        fprintf('  %-10s %9.3g %14.3g %16.3g %9.3g %14.3g %14.3g %14.3g %9.3g\n', ...
            name, w.betaBg, w.z1Bg, w.z2Bg, w.beta, w.ax, w.ay, w.az, priorGamma);
    end
    fprintf('\n');
end

function report = localPrintFieldHealthReport(stageName, fields, sigma, activeMask, config, bounds)
    names = ["EC"; config.names(:)];
    rows = repmat(struct( ...
        'name', "", ...
        'nTotal', 0, ...
        'nActive', 0, ...
        'nComplex', 0, ...
        'nComplexActive', 0, ...
        'nNaN', 0, ...
        'nNaNActive', 0, ...
        'nInf', 0, ...
        'nInfActive', 0, ...
        'nLower', 0, ...
        'nLowerActive', 0, ...
        'minActive', NaN, ...
        'p05Active', NaN, ...
        'medianActive', NaN, ...
        'p95Active', NaN, ...
        'maxActive', NaN), numel(names), 1);

    activeMask = activeMask(:);
    fprintf('\nField numeric health report: %s\n', stageName);
    fprintf('  %-10s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %12s %12s %12s %12s %12s\n', ...
        'field', 'n', 'active', 'complex', 'cxAct', 'NaN', 'nanAct', 'Inf', 'infAct', ...
        'low', 'lowAct', 'minAct', 'p05Act', 'medianAct', 'p95Act', 'maxAct');

    for i = 1:numel(names)
        name = char(names(i));
        if strcmp(name, 'EC')
            vals = sigma(:);
        else
            vals = fields.(name)(:);
        end
        rows(i) = localFieldHealthRow(name, vals, activeMask, bounds);
        fprintf('  %-10s %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %12.4g %12.4g %12.4g %12.4g %12.4g\n', ...
            char(rows(i).name), rows(i).nTotal, rows(i).nActive, rows(i).nComplex, ...
            rows(i).nComplexActive, rows(i).nNaN, rows(i).nNaNActive, rows(i).nInf, rows(i).nInfActive, ...
            rows(i).nLower, rows(i).nLowerActive, ...
            rows(i).minActive, rows(i).p05Active, ...
            rows(i).medianActive, rows(i).p95Active, rows(i).maxActive);
    end
    fprintf('\n');

    report = struct();
    report.stageName = string(stageName);
    report.rows = rows;
end

function localPrintPropertyLayeringReport(stageName, mesh, fields, config)
    nxyz = mesh.gridSize - 1;
    fprintf('\nProperty layering diagnostic: %s\n', stageName);
    fprintf('  %-10s %16s %16s\n', 'property', 'layerMeanRange', 'maxHorizStd');
    for i = 1:numel(config.names)
        name = char(config.names(i));
        vals3 = reshape(fields.(name), nxyz);
        layerMeans = squeeze(mean(mean(vals3, 1), 2));
        horizStd = zeros(nxyz(3), 1);
        for iz = 1:nxyz(3)
            layerVals = vals3(:, :, iz);
            finiteLayerVals = layerVals(isfinite(layerVals));
            if isempty(finiteLayerVals)
                horizStd(iz) = NaN;
            else
                horizStd(iz) = std(finiteLayerVals(:));
            end
        end
        fprintf('  %-10s %16.4g %16.4g\n', ...
            name, max(layerMeans) - min(layerMeans), max(horizStd));
    end
    fprintf('\n');
end

function row = localFieldHealthRow(name, vals, activeMask, bounds)
    vals = vals(:);
    isComplexVal = abs(imag(vals)) > 0;
    valsReal = real(vals);
    activeVals = valsReal(activeMask);
    finiteActive = activeVals(isfinite(activeVals));
    if isempty(finiteActive)
        stats = NaN(1, 5);
    else
        stats = [min(finiteActive), prctile(finiteActive, 5), median(finiteActive), ...
            prctile(finiteActive, 95), max(finiteActive)];
    end

    row = struct();
    row.name = string(name);
    row.nTotal = numel(vals);
    row.nActive = nnz(activeMask);
    row.nComplex = nnz(isComplexVal);
    row.nComplexActive = nnz(isComplexVal(activeMask));
    row.nNaN = nnz(isnan(valsReal));
    row.nNaNActive = nnz(isnan(valsReal(activeMask)));
    row.nInf = nnz(isinf(valsReal));
    row.nInfActive = nnz(isinf(valsReal(activeMask)));
    row.nLower = 0;
    row.nLowerActive = 0;
    if isfield(bounds, name)
        lim = bounds.(name);
        lowerBound = lim(1);
        tol = max(1e-9, 1e-8 * max(1, abs(lowerBound)));
        row.nLower = nnz(valsReal <= lowerBound + tol);
        row.nLowerActive = nnz(valsReal(activeMask) <= lowerBound + tol);
    end
    row.minActive = stats(1);
    row.p05Active = stats(2);
    row.medianActive = stats(3);
    row.p95Active = stats(4);
    row.maxActive = stats(5);
end

function report = localPrintPropertyResidualContributions(mesh, config, fieldsCurrent, refFields, truthSettings, activeMask, weights)
    actualFields = evaluate_soil_true_model_at_points(mesh.elementCenters, refFields, truthSettings);
    activeMask = activeMask(:);
    propNames = config.names(:);
    rows = repmat(struct( ...
        'property', "", ...
        'beta', NaN, ...
        'rmsRaw', NaN, ...
        'rmsRelative', NaN, ...
        'totalRelativeEnergy', NaN, ...
        'weightedRelativeEnergy', NaN, ...
        'suggestedBetaMultiplier', NaN), numel(propNames), 1);

    relEnergy = NaN(numel(propNames), 1);
    for i = 1:numel(propNames)
        name = char(propNames(i));
        currentVals = fieldsCurrent.(name)(activeMask);
        actualVals = actualFields.(name)(activeMask);
        rawResidual = currentVals - actualVals;
        relResidual = rawResidual ./ max(abs(actualVals), 1e-9);

        beta = NaN;
        if isfield(weights, name) && isfield(weights.(name), 'beta')
            beta = weights.(name).beta;
        end
        relEnergy(i) = sum(relResidual.^2);

        rows(i).property = string(name);
        rows(i).beta = beta;
        rows(i).rmsRaw = sqrt(mean(rawResidual.^2));
        rows(i).rmsRelative = sqrt(mean(relResidual.^2));
        rows(i).totalRelativeEnergy = relEnergy(i);
        rows(i).weightedRelativeEnergy = beta * relEnergy(i);
    end

    finiteEnergy = relEnergy(isfinite(relEnergy) & relEnergy > 0);
    if isempty(finiteEnergy)
        referenceEnergy = NaN;
    else
        referenceEnergy = median(finiteEnergy);
    end
    for i = 1:numel(rows)
        rows(i).suggestedBetaMultiplier = rows(i).totalRelativeEnergy / max(referenceEnergy, eps);
    end

    fprintf('\nSoil-property residual contribution report\n');
    fprintf('  Residuals are computed over active soil cells only: (current - synthetic actual) / actual.\n');
    fprintf('  betaMultiplier is relative to the median total relative residual energy; use it as a tuning hint, not a rule.\n');
    fprintf('  %-10s %10s %12s %12s %14s %16s %12s\n', ...
        'property', 'beta', 'rawRMS', 'relRMS', 'relEnergy', 'beta*relEnergy', 'betaMult');
    for i = 1:numel(rows)
        fprintf('  %-10s %10.3g %12.3g %12.3g %14.3g %16.3g %12.3g\n', ...
            char(rows(i).property), rows(i).beta, rows(i).rmsRaw, rows(i).rmsRelative, ...
            rows(i).totalRelativeEnergy, rows(i).weightedRelativeEnergy, rows(i).suggestedBetaMultiplier);
    end
    fprintf('\n');

    report = struct();
    report.activeCellCount = nnz(activeMask);
    report.referenceRelativeEnergy = referenceEnergy;
    report.rows = rows;
end

function report = localPrintPriorFitReport(config, fieldsCurrent, prior, weights)
    report = struct('rows', []);
    if isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'hasMeshMapping') || ~prior.hasMeshMapping
        fprintf('\nPrior fit report: no mapped prior samples available.\n');
        return;
    end
    if ~isfield(weights, 'prior') || isempty(weights.prior)
        fprintf('\nPrior fit report: no prior gamma weights configured.\n');
        return;
    end

    propNames = string(fieldnames(weights.prior));
    rows = repmat(struct( ...
        'property', "", ...
        'gamma', NaN, ...
        'nSamples', 0, ...
        'meanAbsResidual', NaN, ...
        'maxAbsResidual', NaN, ...
        'rmsResidual', NaN), numel(propNames), 1);

    fprintf('\nPrior fit at sampled cells\n');
    fprintf('  %-10s %12s %10s %14s %14s %14s\n', ...
        'property', 'gamma', 'n', 'meanAbs', 'maxAbs', 'rms');
    used = false(numel(propNames), 1);
    for i = 1:numel(propNames)
        name = propNames(i);
        gamma = weights.prior.(char(name));
        if ~(isnumeric(gamma) && isscalar(gamma) && gamma > 0)
            continue;
        end
        if ~ismember(name, config.names) || ~isfield(prior, char(name)) || ~isfield(fieldsCurrent, char(name))
            continue;
        end
        sampleIdx = prior.cellIndex(:);
        currentVals = fieldsCurrent.(char(name))(sampleIdx);
        targetVals = prior.(char(name))(:);
        residual = currentVals - targetVals;

        rows(i).property = name;
        rows(i).gamma = gamma;
        rows(i).nSamples = numel(residual);
        rows(i).meanAbsResidual = mean(abs(residual), 'omitnan');
        rows(i).maxAbsResidual = max(abs(residual), [], 'omitnan');
        rows(i).rmsResidual = sqrt(mean(residual.^2, 'omitnan'));
        used(i) = true;
        fprintf('  %-10s %12.4g %10d %14.4g %14.4g %14.4g\n', ...
            name, gamma, rows(i).nSamples, rows(i).meanAbsResidual, ...
            rows(i).maxAbsResidual, rows(i).rmsResidual);
    end
    if ~any(used)
        fprintf('  No prior properties were active. Check prior gamma values and prior CSV columns.\n');
        report.rows = rows([]);
    else
        report.rows = rows(used);
    end
end

function mOut = localApplyHardPriorSamples(mIn, config, prior, hardPriorSamples)
    mOut = mIn;
    if ~hardPriorSamples
        return;
    end
    if isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'hasMeshMapping') || ~prior.hasMeshMapping || ...
            ~isfield(prior, 'cellIndex')
        return;
    end
    sampleIdx = prior.cellIndex(:);
    for i = 1:numel(config.names)
        name = config.names(i);
        key = char(name);
        if ~isfield(prior, key)
            continue;
        end
        idxAll = config.index.(key);
        targetVals = prior.(key)(:);
        if numel(targetVals) ~= numel(sampleIdx)
            continue;
        end
        mOut(idxAll(sampleIdx)) = soil_property_to_inversion_var(name, targetVals, config);
    end
end

function [mCurrent, sigmaCurrent, fieldsCurrent, iterationLog, summary] = localRunPropertyEcCalibrationProjected( ...
        mesh, sigmaTarget, mStart, mRef, weights, prior, bounds, config, aux, ecModel, projection, activeMask, maxIter, lineSearchSteps, stageName, hardPriorSamples)

    if nargin < 15 || strlength(string(stageName)) == 0
        stageName = "property";
    end
    if nargin < 16
        hardPriorSamples = false;
    end
    activeMask = activeMask(:);
    if numel(activeMask) ~= config.nCells
        error('activeMask must have one entry per inversion cell.');
    end

    mCurrent = localApplyHardPriorSamples(mStart, config, prior, hardPriorSamples);
    [sigmaCurrent, ~, fieldsCurrent] = predict_soil_ec_physics(mCurrent, config, aux, ecModel);
    yTargetAll = log10(max(sigmaTarget, realmin));
    yTarget = yTargetAll(activeMask);
    dataStd = max(0.02, 0.05 * max(abs(yTarget), 1e-3));
    Wsigma = spdiags(1 ./ dataStd, 0, numel(yTarget), numel(yTarget));

    objectiveCurrent = localEvaluatePropertyEcObjective(mCurrent, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaCurrent, activeMask);
    objectiveStart = objectiveCurrent.total;
    iterationLog = repmat(struct('iter',0,'objective',0,'phiData',0,'phiReg',0,'phiTexture',0,'phiPrior',0,'normalizedDataMisfit',0,'stepNorm',0,'stepLength',0,'accepted',false), maxIter, 1);
    exitReason = "maxIter";
    lastStepLength = 0;

    fprintf('\n%s solve\n', localTitleCase(stageName));
    fprintf('  Start: obj=%.4g phiData=%.4g normMisfit=%.4g phiReg=%.4g phiPrior=%.4g phiECprior=%.4g\n', ...
        objectiveCurrent.total, objectiveCurrent.phiData, objectiveCurrent.normalizedDataMisfit, ...
        objectiveCurrent.phiReg, objectiveCurrent.phiPrior, objectiveCurrent.phiPriorEC);

    for iter = 1:maxIter
        fprintf('  Iter %d/%d: building sensitivities and regularization ... ', iter, maxIter);
        [sigmaCurrent, JsoilSigma, fieldsCurrent] = predict_soil_ec_physics(mCurrent, config, aux, ecModel);
        JlogSigma = spdiags(1 ./ (log(10) * max(sigmaCurrent, realmin)), 0, config.nCells, config.nCells) * JsoilSigma;
        [R, rhsReg] = build_property_regularization(config, mesh, mCurrent, mRef, weights, prior);
        [Rec, rhsEc] = build_property_ec_prior_penalty(sigmaCurrent, JsoilSigma, prior, weights);
        R = [R; Rec];
        rhsReg = [rhsReg; rhsEc];

        [Jproj, Rproj, projectionInfo] = localProjectPropertySystem(JlogSigma(activeMask, :), R, projection);
        rData = yTarget - log10(max(sigmaCurrent(activeMask), realmin));
        fprintf('solving %dx%d augmented system ... ', size(Jproj, 2), size(Jproj, 2));
        [deltaProjected, solveInfo] = solve_property_update(Jproj, rData, Wsigma, Rproj, rhsReg);
        deltaM = localExpandPropertyUpdate(deltaProjected, projection);
        fprintf('update solved (mode=%s, rows=%d, cols=%d, damp=%.2g, stepNorm=%.3g)\n', ...
            char(projectionInfo.mode), solveInfo.systemRows, solveInfo.systemCols, ...
            solveInfo.damping, norm(deltaM) / max(norm(mCurrent), eps));

        accepted = false;
        best = struct('objectiveValue', Inf);
        for stepLength = lineSearchSteps
            fprintf('    trial step %.3g ... ', stepLength);
            [mTrial, fieldsTrial, boundInfo] = apply_property_bounds(mCurrent + stepLength * deltaM, config, bounds);
            mTrial = localApplyHardPriorSamples(mTrial, config, prior, hardPriorSamples);
            [mTrial, fieldsTrial, boundInfoPinned] = apply_property_bounds(mTrial, config, bounds);
            boundInfo.nClamped = boundInfo.nClamped + boundInfoPinned.nClamped;
            [sigmaTrial, ~, fieldsTrial] = predict_soil_ec_physics(mTrial, config, aux, ecModel);
            objectiveTrial = localEvaluatePropertyEcObjective(mTrial, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaTrial, activeMask);
            fprintf('obj=%.4g normMisfit=%.4g phiPrior=%.4g phiECprior=%.4g clamped=%d', ...
                objectiveTrial.total, objectiveTrial.normalizedDataMisfit, ...
                objectiveTrial.phiPrior, objectiveTrial.phiPriorEC, boundInfo.nClamped);
            if objectiveTrial.total < best.objectiveValue
                accepted = true;
                best.m = mTrial;
                best.sigma = sigmaTrial;
                best.fields = fieldsTrial;
                best.objective = objectiveTrial;
                best.objectiveValue = objectiveTrial.total;
                best.stepLength = stepLength;
                fprintf(' best-so-far');
            end
            fprintf('\n');
        end

        if accepted && best.objectiveValue >= objectiveCurrent.total
            accepted = false;
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
            fprintf('  Iter %d stopped: no trial step improved the objective.\n', iter);
            break;
        end
        fprintf('    selected best trial step %.3g\n', best.stepLength);

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
        fprintf('  Iter %d accepted: step=%.3g obj=%.4g phiData=%.4g normMisfit=%.4g phiReg=%.4g phiPrior=%.4g phiECprior=%.4g\n', ...
            iter, best.stepLength, objectiveCurrent.total, objectiveCurrent.phiData, ...
            objectiveCurrent.normalizedDataMisfit, objectiveCurrent.phiReg, ...
            objectiveCurrent.phiPrior, objectiveCurrent.phiPriorEC);
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
    summary.priorRows = localCountPriorRows(prior);
    summary.nPriorRowsUsed = summary.priorRows;

    fprintf('  Finished %s: exit=%s iterations=%d obj=%.4g normMisfit=%.4g\n\n', ...
        stageName, summary.exitReason, summary.iterationsCompleted, ...
        summary.objectiveAfter, summary.normalizedDataMisfitAfter);
end

function [Jproj, Rproj, info] = localProjectPropertySystem(J, R, projection)
    mode = string(projection.mode);
    switch mode
        case "matrix"
            P = projection.matrix;
            Jproj = J * P;
            Rproj = R * P;
        case "removeDepthBasis"
            Jproj = J;
            Rproj = R;
        case "removeLayerwiseMean"
            Jproj = J;
            Rproj = R;
        otherwise
            error('Unsupported property projection mode: %s.', mode);
    end
    info = struct('mode', mode);
end

function deltaM = localExpandPropertyUpdate(deltaProjected, projection)
    mode = string(projection.mode);
    switch mode
        case "matrix"
            deltaM = projection.matrix * deltaProjected;
        case "removeDepthBasis"
            deltaM = deltaProjected;
            Q = projection.Q;
            nCells = projection.nCells;
            for i = 1:projection.nProps
                idx = ((i - 1) * nCells + 1):(i * nCells);
                du = deltaM(idx);
                deltaM(idx) = du - Q * (Q' * du);
            end
        case "removeLayerwiseMean"
            deltaM = deltaProjected;
            nCells = projection.nCells;
            nxyz = projection.gridSize;
            for i = 1:projection.nProps
                idx = ((i - 1) * nCells + 1):(i * nCells);
                du3 = reshape(deltaM(idx), nxyz);
                layerMean = mean(mean(du3, 1), 2);
                du3 = du3 - repmat(layerMean, [nxyz(1), nxyz(2), 1]);
                deltaM(idx) = du3(:);
            end
            if isfield(projection, 'updateMask')
                deltaM(~projection.updateMask) = 0;
            end
        otherwise
            error('Unsupported property projection mode: %s.', mode);
    end
end

function [mask, resolvedNames] = localPropertyUpdateMask(config, updateNames)
    names = string(updateNames(:));
    if any(strcmpi(names, "all"))
        mask = true(config.totalSize, 1);
        resolvedNames = config.names(:);
        return;
    end

    mask = false(config.totalSize, 1);
    resolvedNames = strings(0, 1);
    for i = 1:numel(names)
        hit = config.names(strcmpi(config.names, names(i)));
        if isempty(hit)
            error('Unknown propertyAnomalyUpdateNames entry: %s', names(i));
        end
        key = char(hit(1));
        mask(config.index.(key)) = true;
        resolvedNames(end + 1, 1) = hit(1); %#ok<AGROW>
    end
    resolvedNames = unique(resolvedNames, 'stable');
end

function s = localTitleCase(s)
    s = char(string(s));
    if isempty(s)
        return;
    end
    s(1) = upper(s(1));
end

function info = localEvaluatePropertyEcObjective(m, sigmaTarget, Wsigma, config, mesh, mRef, weights, prior, sigmaCurrent, activeMask)
    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    fields = unpack_soil_fields(m, config);
    if nargin < 10 || isempty(activeMask)
        activeMask = true(config.nCells, 1);
    end
    activeMask = activeMask(:);
    dataResidual = Wsigma * (log10(max(sigmaTarget(activeMask), realmin)) - log10(max(sigmaCurrent(activeMask), realmin)));
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
        duRef = u - uRef;
        w = weights.(char(name));
        phiReg = phiReg + w.beta * sum(duRef.^2);
        phiReg = phiReg + w.ax * sum((Dx * duRef).^2) + w.ay * sum((Dy * duRef).^2) + w.az * sum((Dz * duRef).^2);
        phiReg = phiReg + w.cx * sum((Dxx * duRef).^2) + w.cy * sum((Dyy * duRef).^2) + w.cz * sum((Dzz * duRef).^2);
        if isfield(w, 'depthBasisGamma') && w.depthBasisGamma > 0
            phiDepthBasis = phiDepthBasis + w.depthBasisGamma * sum((Pperp * duRef).^2);
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
    info.normalizedDataMisfit = phiData / max(1, nnz(activeMask));
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

function aux = localMakeAux(nCells)
    aux = struct();
    aux.pH_CaCl2 = 6.4 * ones(nCells, 1);
    aux.pH_H2O = 7.0 * ones(nCells, 1);
    aux.CaCO3 = 2.0 * ones(nCells, 1);
    aux.Coarse = 5.0 * ones(nCells, 1);
end

function n = localCountPriorRows(prior)
    if isempty(prior) || ~isfield(prior, 'nSamples')
        n = 0;
    else
        n = prior.nSamples;
    end
end
