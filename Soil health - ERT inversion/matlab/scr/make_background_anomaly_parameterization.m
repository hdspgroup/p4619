function param = make_background_anomaly_parameterization(mesh, options)
%MAKE_BACKGROUND_ANOMALY_PARAMETERIZATION Build u = B*u_bg + u_anom mapping.

    arguments
        mesh struct
        options.BackgroundMode (1,1) string = "layerwise" % constant | layerwise | gaussian3
    end

    n = mesh.gridSize - 1;
    nCells = prod(n);
    nZ = n(3);
    zCenters = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));

    idx = reshape(1:nCells, n);

    switch lower(options.BackgroundMode)
        case "constant"
            basisZ = ones(nZ, 1);

        case "layerwise"
            basisZ = speye(nZ);

        case "gaussian3"
            basisZ = localGaussianBasis(zCenters, 3);

        otherwise
            error('Unsupported BackgroundMode: %s', options.BackgroundMode);
    end

    nBg = size(basisZ, 2);
    B = localLiftBasisToCells(idx, basisZ, nCells);

    % Depthwise averaging operator for the anomaly: one mean value per z
    % level, used to keep the anomaly as a true perturbation around the
    % layered background.
    avgRows = [];
    avgCols = [];
    avgVals = [];
    for iz = 1:nZ
        ids = idx(:, :, iz);
        ids = ids(:);
        nHere = numel(ids);
        avgRows = [avgRows; iz * ones(nHere, 1)]; %#ok<AGROW>
        avgCols = [avgCols; ids]; %#ok<AGROW>
        avgVals = [avgVals; (1 / nHere) * ones(nHere, 1)]; %#ok<AGROW>
    end
    Amean = sparse(avgRows, avgCols, avgVals, nZ, nCells);

    G = [B, speye(nCells)];

    param = struct();
    param.nCells = nCells;
    param.nZ = nZ;
    param.nBg = nBg;
    param.backgroundMode = string(options.BackgroundMode);
    param.zCenters = zCenters(:);
    param.basisZ = basisZ;
    param.B = B;
    param.Amean = Amean;
    param.G = G;
    param.idxBg = 1:nBg;
    param.idxAnom = nBg + (1:nCells);
end

function B = localLiftBasisToCells(idx, basisZ, nCells)
    nZ = size(basisZ, 1);
    nBg = size(basisZ, 2);

    rowIdx = [];
    colIdx = [];
    vals = [];
    for iz = 1:nZ
        ids = idx(:, :, iz);
        ids = ids(:);
        nHere = numel(ids);
        for ib = 1:nBg
            rowIdx = [rowIdx; ids]; %#ok<AGROW>
            colIdx = [colIdx; ib * ones(nHere, 1)]; %#ok<AGROW>
            vals = [vals; basisZ(iz, ib) * ones(nHere, 1)]; %#ok<AGROW>
        end
    end
    B = sparse(rowIdx, colIdx, vals, nCells, nBg);
end

function basisZ = localGaussianBasis(zCenters, nBasis)
    zMin = min(zCenters);
    zMax = max(zCenters);
    centers = linspace(zMin, zMax, nBasis);
    spacing = max((zMax - zMin) / max(nBasis - 1, 1), eps);
    sigma = 1.25 * spacing;

    basisZ = zeros(numel(zCenters), nBasis);
    for i = 1:nBasis
        basisZ(:, i) = exp(-0.5 * ((zCenters - centers(i)) / sigma) .^ 2);
    end

    rowSums = sum(basisZ, 2);
    rowSums(rowSums == 0) = 1;
    basisZ = basisZ ./ rowSums;
end
