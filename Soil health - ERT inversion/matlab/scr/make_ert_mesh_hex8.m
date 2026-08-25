function mesh = make_ert_mesh_hex8(options)
%MAKE_ERT_MESH_HEX8 Create a structured rectilinear mesh for ERT tests.
%
% Historically this returned a hex8 FEM mesh. It now also carries the
% cell-centered geometry needed by the finite-volume ERT solver while
% preserving the familiar structured-mesh fields used throughout the
% project.

    arguments
        options.XLim (1,2) double = [-45, 45]
        options.YLim (1,2) double = [-45, 45]
        options.ZLim (1,2) double = [0, 45]
        options.Spacing (1,1) double = 2.5
        options.SpacingXYZ (1,3) double = [NaN, NaN, NaN]
        options.XCoords double = []
        options.YCoords double = []
        options.ZCoords double = []
    end

    if ~isempty(options.XCoords)
        x = localValidateAxis(options.XCoords, options.XLim, 'XCoords');
    else
        if all(isfinite(options.SpacingXYZ))
            dx = options.SpacingXYZ(1);
        else
            dx = options.Spacing;
        end
        x = options.XLim(1):dx:options.XLim(2);
        if x(end) ~= options.XLim(2)
            x = [x, options.XLim(2)];
        end
    end

    if ~isempty(options.YCoords)
        y = localValidateAxis(options.YCoords, options.YLim, 'YCoords');
    else
        if all(isfinite(options.SpacingXYZ))
            dy = options.SpacingXYZ(2);
        else
            dy = options.Spacing;
        end
        y = options.YLim(1):dy:options.YLim(2);
        if y(end) ~= options.YLim(2)
            y = [y, options.YLim(2)];
        end
    end

    if ~isempty(options.ZCoords)
        z = localValidateAxis(options.ZCoords, options.ZLim, 'ZCoords');
    else
        if all(isfinite(options.SpacingXYZ))
            dz = options.SpacingXYZ(3);
        else
            dz = options.Spacing;
        end
        z = options.ZLim(1):dz:options.ZLim(2);
        if z(end) ~= options.ZLim(2)
            z = [z, options.ZLim(2)];
        end
    end

    [X, Y, Z] = ndgrid(x, y, z);
    nodes = [X(:), Y(:), Z(:)];

    nx = numel(x);
    ny = numel(y);
    nz = numel(z);
    nElements = (nx - 1) * (ny - 1) * (nz - 1);
    elements = zeros(nElements, 8);

    e = 0;
    for k = 1:nz-1
        for j = 1:ny-1
            for i = 1:nx-1
                e = e + 1;
                elements(e, :) = [ ...
                    sub2ind([nx, ny, nz], i,   j,   k), ...
                    sub2ind([nx, ny, nz], i+1, j,   k), ...
                    sub2ind([nx, ny, nz], i+1, j+1, k), ...
                    sub2ind([nx, ny, nz], i,   j+1, k), ...
                    sub2ind([nx, ny, nz], i,   j,   k+1), ...
                    sub2ind([nx, ny, nz], i+1, j,   k+1), ...
                    sub2ind([nx, ny, nz], i+1, j+1, k+1), ...
                    sub2ind([nx, ny, nz], i,   j+1, k+1)];
            end
        end
    end

    elementCenters = zeros(nElements, 3);
    for e = 1:nElements
        elementCenters(e, :) = mean(nodes(elements(e, :), :), 1);
    end

    xc = 0.5 * (x(1:end-1) + x(2:end));
    yc = 0.5 * (y(1:end-1) + y(2:end));
    zc = 0.5 * (z(1:end-1) + z(2:end));
    [Xc, Yc, Zc] = ndgrid(xc, yc, zc);

    mesh = struct();
    mesh.elementType = "fvm_rect";
    mesh.discretization = "fvm";
    mesh.meshStyle = "rect_refined";
    mesh.nodes = nodes;
    mesh.elements = elements;
    mesh.elementCenters = [Xc(:), Yc(:), Zc(:)];
    mesh.x = x;
    mesh.y = y;
    mesh.z = z;
    mesh.xCenters = xc(:);
    mesh.yCenters = yc(:);
    mesh.zCenters = zc(:);
    mesh.dx = diff(x(:)');
    mesh.dy = diff(y(:)');
    mesh.dz = diff(z(:)');
    mesh.gridSize = [nx, ny, nz];
    mesh.nCellsXYZ = [nx - 1, ny - 1, nz - 1];
    mesh.nNodes = size(nodes, 1);
    mesh.nElements = nElements;
    mesh.nCells = nElements;
    [DX, DY, DZ] = ndgrid(mesh.dx, mesh.dy, mesh.dz);
    mesh.cellVolumes = (DX .* DY .* DZ);
    mesh.cellVolumes = mesh.cellVolumes(:);
end

function axisVals = localValidateAxis(axisVals, axisLim, argName)
%LOCALVALIDATEAXIS Validate an explicit structured-mesh axis.

    axisVals = axisVals(:)';
    if numel(axisVals) < 2
        error('%s must contain at least two coordinates.', argName);
    end
    if any(diff(axisVals) <= 0)
        error('%s must be strictly increasing.', argName);
    end
    if abs(axisVals(1) - axisLim(1)) > 1e-10 || abs(axisVals(end) - axisLim(2)) > 1e-10
        error('%s must start at %.12g and end at %.12g.', argName, axisLim(1), axisLim(2));
    end
end
