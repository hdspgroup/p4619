function result = solve_ert_forward(mesh, survey, A, options)
%SOLVE_ERT_FORWARD Solve a cell-centered finite-volume ERT forward problem.
%
% Unknowns are cell-centered potentials on a rectilinear mesh. Electrodes
% are tied to the nearest shallow cell center and the measurement voltages
% are extracted from those cell potentials.

    arguments
        mesh struct
        survey struct
        A double
        options.ReferenceNode (1,1) double = NaN
    end

    if ~isfield(mesh, 'xCenters') || ~isfield(mesh, 'yCenters') || ~isfield(mesh, 'zCenters')
        error('FVM forward solve requires xCenters, yCenters, and zCenters in the mesh.');
    end

    nx = numel(mesh.xCenters);
    ny = numel(mesh.yCenters);
    nz = numel(mesh.zCenters);
    nCells = nx * ny * nz;
    if ~isequal(size(A), [nCells, nCells])
        error('System matrix must be nCells x nCells for the provided mesh.');
    end

    elecCells = localMapElectrodesToTopCells(mesh, survey.electrodes);

    sourcePairs = unique(survey.quads(:, 1:2), 'rows', 'stable');
    nSources = size(sourcePairs, 1);
    sourcePairIndex = zeros(survey.nMeasurements, 1);
    for i = 1:survey.nMeasurements
        sourcePairIndex(i) = find(sourcePairs(:, 1) == survey.quads(i, 1) & sourcePairs(:, 2) == survey.quads(i, 2), 1, 'first');
    end

    rhs = sparse(nCells, nSources);
    cellVolumes = mesh.cellVolumes(:);
    for i = 1:nSources
        aCell = elecCells(sourcePairs(i, 1));
        bCell = elecCells(sourcePairs(i, 2));
        rhs(aCell, i) = rhs(aCell, i) + survey.current / max(cellVolumes(aCell), realmin);
        rhs(bCell, i) = rhs(bCell, i) - survey.current / max(cellVolumes(bCell), realmin);
    end

    if isnan(options.ReferenceNode)
        refCell = nCells;
    else
        refCell = max(1, min(nCells, round(options.ReferenceNode)));
    end

    Asolve = A;
    Asolve(refCell, :) = 0;
    Asolve(:, refCell) = 0;
    Asolve(refCell, refCell) = 1;
    rhs(refCell, :) = 0;

    phiBySource = Asolve \ rhs;

    voltages = zeros(survey.nMeasurements, 1);
    mCells = elecCells(survey.quads(:, 3));
    nCellsIdx = elecCells(survey.quads(:, 4));
    for i = 1:survey.nMeasurements
        src = sourcePairIndex(i);
        voltages(i) = phiBySource(mCells(i), src) - phiBySource(nCellsIdx(i), src);
    end

    result = struct();
    result.voltage = voltages;
    result.voltages = voltages;
    result.phiBySource = phiBySource;
    result.sourcePairs = sourcePairs;
    result.sourcePairIndex = sourcePairIndex;
    result.electrodeCellIndex = elecCells;
    result.referenceNode = refCell;
end

function elecCells = localMapElectrodesToTopCells(mesh, electrodes)
    xc = mesh.xCenters(:)';
    yc = mesh.yCenters(:)';
    topK = 1;
    elecCells = zeros(size(electrodes, 1), 1);
    for i = 1:size(electrodes, 1)
        [~, ix] = min(abs(xc - electrodes(i, 1)));
        [~, iy] = min(abs(yc - electrodes(i, 2)));
        elecCells(i) = sub2ind([numel(xc), numel(yc), numel(mesh.zCenters)], ix, iy, topK);
    end
end
