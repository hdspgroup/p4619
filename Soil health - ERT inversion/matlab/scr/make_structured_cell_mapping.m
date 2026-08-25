function [fineFromCoarse, coarseFromFine] = make_structured_cell_mapping(coarseMesh, fineMesh)
%MAKE_STRUCTURED_CELL_MAPPING Map coarse cell values onto a fine structured mesh.
%
% fineFromCoarse maps coarse cell values to fine cells by containment.
% coarseFromFine averages fine cell values back to the coarse grid.

    xBin = discretize(fineMesh.elementCenters(:, 1), coarseMesh.x);
    yBin = discretize(fineMesh.elementCenters(:, 2), coarseMesh.y);
    zBin = discretize(fineMesh.elementCenters(:, 3), coarseMesh.z);

    if any(isnan(xBin) | isnan(yBin) | isnan(zBin))
        error('Fine mesh contains element centers outside the coarse mesh bounds.');
    end

    coarseGrid = coarseMesh.gridSize - 1;
    coarseIdx = sub2ind(coarseGrid, xBin, yBin, zBin);

    nFine = fineMesh.nElements;
    nCoarse = coarseMesh.nElements;

    fineRows = (1:nFine)';
    fineFromCoarse = sparse(fineRows, coarseIdx, 1, nFine, nCoarse);

    counts = full(sum(fineFromCoarse, 1))';
    counts(counts == 0) = 1;
    coarseFromFine = spdiags(1 ./ counts, 0, nCoarse, nCoarse) * fineFromCoarse';
end
