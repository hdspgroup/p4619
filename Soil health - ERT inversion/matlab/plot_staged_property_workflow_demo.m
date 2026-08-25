clear; clc; close all;

% Workflow-style plot for the staged soil-property inversion:
% true properties -> soil-to-EC forward model -> recovered properties.

saveArg = true;
ySlice = 0;
plotProperties = ["N", "OC", "Clay"];

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(scriptDir, 'scr')));

outDir = fullfile(rootDir, 'outputs', 'ert_staged_ec_property_inversion_demo');
propertyFile = fullfile(outDir, 'staged_property_inversion_output.mat');
ecFile = fullfile(outDir, 'staged_ec_inversion_output.mat');
if ~isfile(propertyFile)
    error('Missing staged property output: %s', propertyFile);
end
if ~isfile(ecFile)
    error('Missing staged EC output: %s', ecFile);
end

P = load(propertyFile);
E = load(ecFile);

mesh = P.mesh;
config = P.config;
fieldsCurrent = P.fieldsProp;
ecModel = load_forward_soil_ec_physics_model( ...
    fullfile(rootDir, 'soil_health_ec_models', 'forward_soil_ec_forward_model_coefficients.csv'));

[~, refFields] = make_soil_reference_model(config);
if isfield(E, 'truthSettings')
    truthSettings = E.truthSettings;
else
    truthSettings = [];
end

% Plot the actuals on the same cell mesh used by the property inversion.
% A dense analytic display is attractive, but it can imply resolution that
% the inversion did not actually use.
fieldsTrue = evaluate_soil_true_model_at_points(mesh.elementCenters, refFields, truthSettings);
if isfield(P, 'sigmaTarget') && numel(P.sigmaTarget) == mesh.nElements
    sigmaTrue = P.sigmaTarget(:);
else
    sigmaTrue = localPredictSigmaFromFields(fieldsTrue, config, P, ecModel);
end

fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.06, 0.08, 0.88, 0.78]);
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Soil-property inversion workflow: actual properties -> EC forward model -> recovered properties');

axActual = gobjects(3, 1);
axInverse = gobjects(3, 1);
for i = 1:3
    prop = char(plotProperties(i));
    actualGrid = localExtractXZSlice(mesh, fieldsTrue.(prop), ySlice);
    currentGrid = localExtractXZSlice(mesh, fieldsCurrent.(prop), ySlice);
    climVals = localSharedClim(actualGrid, currentGrid);

    axActual(i) = nexttile((i - 1) * 3 + 1);
    plot_rectilinear_cells(mesh.x, mesh.z, actualGrid);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('%s actual', prop)); xlabel('x (m)'); ylabel('z (m)');

    if i ~= 2
        nexttile((i - 1) * 3 + 2);
        axis off;
    end

    axInverse(i) = nexttile((i - 1) * 3 + 3);
    plot_rectilinear_cells(mesh.x, mesh.z, currentGrid);
    set(gca, 'YDir', 'reverse'); axis tight;
    colormap(gca, jet); clim(climVals); colorbar;
    title(sprintf('%s recovered', prop)); xlabel('x (m)'); ylabel('z (m)');
end

axEc = nexttile(5);
cla(axEc);
set(axEc, 'Visible', 'on');
sigmaTrueXZ = localExtractXZSlice(mesh, sigmaTrue, ySlice);
plot_rectilinear_cells(mesh.x, mesh.z, sigmaTrueXZ);
set(gca, 'YDir', 'reverse'); axis tight;
colormap(gca, jet); colorbar;
title('Forward soil-to-EC model'); xlabel('x (m)'); ylabel('z (m)');

drawnow;
localAddWorkflowArrows(fig, axActual, axEc, axInverse);

if saveArg
    outPng = fullfile(outDir, 'ert_staged_property_workflow_demo.png');
    exportgraphics(fig, outPng, 'Resolution', 220);
    fprintf('Workflow figure written to %s\n', outPng);
end

function sigma = localPredictSigmaFromFields(fields, config, source, ecModel)
    nCells = numel(fields.N);
    cfg = config;
    cfg.nCells = nCells;
    cfg.totalSize = cfg.nProps * nCells;
    index = struct();
    for i = 1:cfg.nProps
        first = (i - 1) * nCells + 1;
        last = i * nCells;
        index.(char(cfg.names(i))) = first:last;
    end
    cfg.index = index;

    aux = localMakeAuxFromSource(source, nCells);
    m = pack_soil_fields(fields, cfg);
    sigma = predict_soil_ec_physics(m, cfg, aux, ecModel);
end

function aux = localMakeAuxFromSource(source, nCells)
    aux = struct();
    defaultAux = struct('pH_CaCl2', 6.4, 'pH_H2O', 7.0, 'CaCO3', 2.0, 'Coarse', 5.0, 'CEC', 10.0);
    names = fieldnames(defaultAux);
    if isfield(source, 'aux')
        sourceAux = source.aux;
    else
        sourceAux = struct();
    end
    for i = 1:numel(names)
        name = names{i};
        if isfield(sourceAux, name)
            val = sourceAux.(name);
            if isscalar(val)
                aux.(name) = repmat(val, nCells, 1);
            else
                aux.(name) = val(1) * ones(nCells, 1);
            end
        else
            aux.(name) = defaultAux.(name) * ones(nCells, 1);
        end
    end
end

function vals2 = localExtractXZSlice(mesh, values, yTarget)
    vals3 = reshape(values, mesh.gridSize - 1);
    yCenters = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    [~, iy] = min(abs(yCenters - yTarget));
    vals2 = squeeze(vals3(:, iy, :))';
end

function climVals = localSharedClim(actualVals, currentVals)
    vals = [actualVals(:); currentVals(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        climVals = [0, 1];
        return;
    end
    climVals = [min(vals), max(vals)];
    if climVals(1) == climVals(2)
        pad = max(abs(climVals(1)), 1) * 0.05;
        climVals = climVals + [-pad, pad];
    end
end

function localAddWorkflowArrows(fig, axActual, axEc, axInverse)
    ecPos = getpixelposition(axEc, true);
    figPos = getpixelposition(fig, true);
    ecLeft = [ecPos(1), ecPos(2) + 0.50 * ecPos(4)] ./ figPos(3:4);
    ecRight = [ecPos(1) + ecPos(3), ecPos(2) + 0.50 * ecPos(4)] ./ figPos(3:4);

    for i = 1:numel(axActual)
        aPos = getpixelposition(axActual(i), true);
        invPos = getpixelposition(axInverse(i), true);

        startLeft = [aPos(1) + aPos(3), aPos(2) + 0.50 * aPos(4)] ./ figPos(3:4);
        stopLeft = ecLeft;
        annotation(fig, 'arrow', [startLeft(1), stopLeft(1)], [startLeft(2), stopLeft(2)], ...
            'LineWidth', 1.4, 'Color', [0.20 0.20 0.20]);

        startRight = ecRight;
        stopRight = [invPos(1), invPos(2) + 0.50 * invPos(4)] ./ figPos(3:4);
        annotation(fig, 'arrow', [startRight(1), stopRight(1)], [startRight(2), stopRight(2)], ...
            'LineWidth', 1.4, 'Color', [0.20 0.20 0.20]);
    end
end
