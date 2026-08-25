function mesh = make_ert_mesh_fvm_octree_style(options)
%MAKE_ERT_MESH_FVM_OCTREE_STYLE Build a dyadically graded rectilinear FVM mesh.
%
% This is "OcTree-style" rather than a full hanging-node octree. Cell sizes
% grow by powers of two away from the refined core, which gives us the same
% practical behavior for the current finite-volume solver while staying on
% an orthogonal rectilinear grid.

    arguments
        options.Electrodes (:,3) double = zeros(0, 3)
        options.CenterXY (1,2) double = [0, 0]
        options.Radius (1,1) double = 40
        options.CoreHalfWidth (1,1) double = 15
        options.MinCellSizeXY (1,1) double = 2.5
        options.MaxCellSizeXY (1,1) double = 20
        options.SurfaceRefineDepth (1,1) double = 3
        options.MinCellSizeZ (1,1) double = 0.5
        options.MaxCellSizeZ (1,1) double = 4
        options.ZMax (1,1) double = 10
        options.OuterPadding (1,1) double = 0
    end

    if ~isempty(options.Electrodes)
        xy = options.Electrodes(:, 1:2);
        centerXY = mean(xy, 1);
        surveyRadius = max(vecnorm(xy - centerXY, 2, 2));
        radius = max(options.Radius, surveyRadius + options.OuterPadding);
        coreHalfWidth = max(options.CoreHalfWidth, 0.6 * surveyRadius);
    else
        centerXY = options.CenterXY;
        radius = options.Radius;
        coreHalfWidth = options.CoreHalfWidth;
    end

    x = centerXY(1) + localBuildDyadicAxis(radius, coreHalfWidth, options.MinCellSizeXY, options.MaxCellSizeXY);
    y = centerXY(2) + localBuildDyadicAxis(radius, coreHalfWidth, options.MinCellSizeXY, options.MaxCellSizeXY);
    z = localBuildDyadicDepthAxis(options.ZMax, options.SurfaceRefineDepth, options.MinCellSizeZ, options.MaxCellSizeZ);

    mesh = make_ert_mesh_hex8( ...
        'XLim', [x(1), x(end)], ...
        'YLim', [y(1), y(end)], ...
        'ZLim', [z(1), z(end)], ...
        'XCoords', x, ...
        'YCoords', y, ...
        'ZCoords', z);
    mesh.meshStyle = "octree_style";
end

function axisEdges = localBuildDyadicAxis(radius, coreHalfWidth, minCell, maxCell)
    if radius <= 0
        axisEdges = [-minCell, 0, minCell];
        return;
    end

    coreHalfWidth = min(max(coreHalfWidth, minCell), radius);
    posEdges = 0;
    current = 0;
    step = minCell;
    while current < radius - 1e-10
        if current >= coreHalfWidth - 1e-10
            step = min(2 * step, maxCell);
        end
        nextEdge = min(current + step, radius);
        posEdges(end + 1) = nextEdge; %#ok<AGROW>
        current = nextEdge;
    end

    axisEdges = unique([-fliplr(posEdges(2:end)), posEdges], 'stable');
end

function zEdges = localBuildDyadicDepthAxis(zMax, refineDepth, minCell, maxCell)
    zMax = max(zMax, minCell);
    refineDepth = min(max(refineDepth, minCell), zMax);

    zEdges = 0;
    current = 0;
    step = minCell;
    while current < zMax - 1e-10
        if current >= refineDepth - 1e-10
            step = min(2 * step, maxCell);
        end
        nextEdge = min(current + step, zMax);
        zEdges(end + 1) = nextEdge; %#ok<AGROW>
        current = nextEdge;
    end
end
