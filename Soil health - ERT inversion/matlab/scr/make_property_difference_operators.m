function [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh)
%MAKE_PROPERTY_DIFFERENCE_OPERATORS Metric-aware first/second derivatives.
%
% First derivatives have units of model-units / m. Second derivatives have
% units of model-units / m^2. This keeps regularization behavior more stable
% when the mesh resolution changes.

    n = mesh.gridSize - 1;
    idx = reshape(1:prod(n), n);

    nCells = prod(n);
    xc = 0.5 * (mesh.x(1:end-1) + mesh.x(2:end));
    yc = 0.5 * (mesh.y(1:end-1) + mesh.y(2:end));
    zc = 0.5 * (mesh.z(1:end-1) + mesh.z(2:end));

    hx = max(diff(xc), eps);
    hy = max(diff(yc), eps);
    hz = max(diff(zc), eps);

    [Hx, ~, ~] = ndgrid(hx, ones(1, n(2)), ones(1, n(3)));
    [~, Hy, ~] = ndgrid(ones(1, n(1)), hy, ones(1, n(3)));
    [~, ~, Hz] = ndgrid(ones(1, n(1)), ones(1, n(2)), hz);

    Dx = localDiffMatrix(idx(1:end-1,:,:)  , idx(2:end,:,:), Hx, nCells);
    Dy = localDiffMatrix(idx(:,1:end-1,:)  , idx(:,2:end,:), Hy, nCells);
    Dz = localDiffMatrix(idx(:,:,1:end-1)  , idx(:,:,2:end), Hz, nCells);

    Dxx = localSecondDiffMatrix(idx(1:end-2,:,:)  , idx(2:end-1,:,:), idx(3:end,:,:), hx(1:end-1), hx(2:end), "x", n, nCells);
    Dyy = localSecondDiffMatrix(idx(:,1:end-2,:)  , idx(:,2:end-1,:), idx(:,3:end,:), hy(1:end-1), hy(2:end), "y", n, nCells);
    Dzz = localSecondDiffMatrix(idx(:,:,1:end-2)  , idx(:,:,2:end-1), idx(:,:,3:end), hz(1:end-1), hz(2:end), "z", n, nCells);
end

function D = localDiffMatrix(leftIdx, rightIdx, h, nCells)
    leftIdx = leftIdx(:);
    rightIdx = rightIdx(:);
    h = max(h(:), eps);
    nRows = numel(leftIdx);
    rr = repelem((1:nRows)', 2, 1);
    cc = [leftIdx; rightIdx];
    vv = [-1 ./ h; 1 ./ h];
    D = sparse(rr, cc, vv, nRows, nCells);
end

function D = localSecondDiffMatrix(leftIdx, midIdx, rightIdx, hLeftAxis, hRightAxis, axisName, n, nCells)
    leftIdx = leftIdx(:);
    midIdx = midIdx(:);
    rightIdx = rightIdx(:);
    nRows = numel(leftIdx);
    if nRows == 0
        D = sparse(0, nCells);
        return;
    end
    hLeft = localExpandSecondSpacing(hLeftAxis, axisName, n);
    hRight = localExpandSecondSpacing(hRightAxis, axisName, n);
    hLeft = max(hLeft(:), eps);
    hRight = max(hRight(:), eps);
    denom = hLeft + hRight;
    vLeft = 2 ./ (hLeft .* denom);
    vMid = -2 .* (1 ./ hLeft + 1 ./ hRight) ./ denom;
    vRight = 2 ./ (hRight .* denom);

    rr = repelem((1:nRows)', 3, 1);
    cc = [leftIdx; midIdx; rightIdx];
    vv = [vLeft; vMid; vRight];
    D = sparse(rr, cc, vv, nRows, nCells);
end

function h = localExpandSecondSpacing(hAxis, axisName, n)
    switch axisName
        case "x"
            [h, ~, ~] = ndgrid(hAxis, ones(1, n(2)), ones(1, n(3)));
        case "y"
            [~, h, ~] = ndgrid(ones(1, n(1)), hAxis, ones(1, n(3)));
        case "z"
            [~, ~, h] = ndgrid(ones(1, n(1)), ones(1, n(2)), hAxis);
        otherwise
            error('Unsupported axisName: %s', axisName);
    end
end
