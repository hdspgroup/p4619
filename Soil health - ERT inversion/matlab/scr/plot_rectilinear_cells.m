function h = plot_rectilinear_cells(xEdges, yEdges, cellValues, varargin)
%PLOT_RECTILINEAR_CELLS Plot cell-centered values on a rectilinear mesh.
%
% xEdges    : 1 x (nx+1) cell-edge coordinates
% yEdges    : 1 x (ny+1) cell-edge coordinates
% cellValues: ny x nx matrix, one value per cell
%
% This uses flat-colored patches, so:
% - axes reflect the true, nonuniform cell geometry
% - each cell keeps its own color without interpolation
%
% Name-value options:
%   'AlphaData' : ny x nx matrix of per-cell opacities in [0, 1]

    parser = inputParser;
    addParameter(parser, 'AlphaData', []);
    parse(parser, varargin{:});
    alphaData = parser.Results.AlphaData;

    if numel(xEdges) < 2 || numel(yEdges) < 2
        error('xEdges and yEdges must contain at least two coordinates.');
    end

    [ny, nx] = size(cellValues);
    if nx ~= numel(xEdges) - 1 || ny ~= numel(yEdges) - 1
        error('cellValues must have size [numel(yEdges)-1, numel(xEdges)-1].');
    end
    if ~isempty(alphaData) && ~isequal(size(alphaData), size(cellValues))
        error('AlphaData must match the size of cellValues.');
    end

    [Xv, Yv] = meshgrid(xEdges, yEdges);
    vertices = [Xv(:), Yv(:)];

    nVertY = numel(yEdges);
    nCells = nx * ny;
    faces = zeros(nCells, 4);
    cdata = zeros(nCells, 1);
    adata = zeros(nCells, 1);

    k = 0;
    for iy = 1:ny
        for ix = 1:nx
            k = k + 1;
            v1 = sub2ind([nVertY, numel(xEdges)], iy, ix);
            v2 = sub2ind([nVertY, numel(xEdges)], iy, ix + 1);
            v3 = sub2ind([nVertY, numel(xEdges)], iy + 1, ix + 1);
            v4 = sub2ind([nVertY, numel(xEdges)], iy + 1, ix);
            faces(k, :) = [v1, v2, v3, v4];
            cdata(k) = cellValues(iy, ix);
            if ~isempty(alphaData)
                adata(k) = alphaData(iy, ix);
            end
        end
    end

    patchArgs = {'Faces', faces, ...
        'Vertices', vertices, ...
        'FaceVertexCData', cdata, ...
        'FaceColor', 'flat', ...
        'EdgeColor', 'none'};

    if ~isempty(alphaData)
        patchArgs = [patchArgs, {'FaceVertexAlphaData', adata, ...
            'FaceAlpha', 'flat', ...
            'AlphaDataMapping', 'none'}];
    end

    h = patch(patchArgs{:});
end
