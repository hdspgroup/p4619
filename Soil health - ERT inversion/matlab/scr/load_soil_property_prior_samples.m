function prior = load_soil_property_prior_samples(filePath, varargin)
%LOAD_SOIL_PROPERTY_PRIOR_SAMPLES Load shallow soil-property prior samples from CSV.
%
% Expected CSV columns:
%   sample_id,x_m,y_m,pH_CaCl2,pH_H2O,point_ec_s_m,N,P,K,OC, ...
%   moisture_pct,clay_pct,silt_pct,sand_pct,CaCO3,Coarse
%
% Example:
%   prior = load_soil_property_prior_samples(csvPath);
%   prior = load_soil_property_prior_samples(csvPath, 'Mesh', mesh);
%
% If a mesh is supplied, each sample is mapped to the nearest mesh cell
% center using a shallow point [x, y, 0], which is a reasonable default for
% shallow soil-core / grab-sample priors.

    parser = inputParser;
    addRequired(parser, 'filePath', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'Mesh', [], @(x) isempty(x) || isstruct(x));
    addParameter(parser, 'TextureSumTolerance', 2.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(parser, filePath, varargin{:});

    filePath = string(parser.Results.filePath);
    mesh = parser.Results.Mesh;
    textureTol = parser.Results.TextureSumTolerance;

    if ~isfile(filePath)
        error('Prior sample file not found: %s', filePath);
    end

    tbl = readtable(filePath, 'TextType', 'string');
    expected = ["sample_id","x_m","y_m","pH_CaCl2","pH_H2O","point_ec_s_m","N","P","K","OC","moisture_pct","clay_pct","silt_pct","sand_pct","CaCO3","Coarse"];
    missing = expected(~ismember(expected, string(tbl.Properties.VariableNames)));
    if ~isempty(missing)
        error('Prior sample file is missing required columns: %s', strjoin(missing, ', '));
    end

    tbl = tbl(:, cellstr(expected));
    tbl.sample_id = string(tbl.sample_id);

    numericCols = ["x_m","y_m","pH_CaCl2","pH_H2O","point_ec_s_m","N","P","K","OC","moisture_pct","clay_pct","silt_pct","sand_pct","CaCO3","Coarse"];
    for name = numericCols
        vals = tbl.(name);
        if ~isnumeric(vals)
            error('Column %s must be numeric.', name);
        end
        if any(~isfinite(vals))
            error('Column %s contains NaN or Inf values.', name);
        end
    end

    textureSum = tbl.clay_pct + tbl.silt_pct + tbl.sand_pct;
    if any(abs(textureSum - 100) > textureTol)
        warning('Some prior samples have clay+silt+sand not close to 100 within tolerance %.2f.', textureTol);
    end

    n = height(tbl);
    prior = struct();
    prior.filePath = char(filePath);
    prior.table = tbl;
    prior.nSamples = n;
    prior.sampleId = tbl.sample_id;
    prior.xy = [tbl.x_m, tbl.y_m];
    prior.pH = tbl.pH_CaCl2;
    prior.pH_CaCl2 = tbl.pH_CaCl2;
    prior.pH_H2O = tbl.pH_H2O;
    prior.PointEC = tbl.point_ec_s_m;
    prior.N = tbl.N;
    prior.P = tbl.P;
    prior.K = tbl.K;
    prior.OC = tbl.OC;
    prior.Moisture = tbl.moisture_pct;
    prior.Clay = tbl.clay_pct;
    prior.Silt = tbl.silt_pct;
    prior.Sand = tbl.sand_pct;
    prior.CaCO3 = tbl.CaCO3;
    prior.Coarse = tbl.Coarse;
    prior.textureSum = textureSum;
    prior.propertyNames = ["pH_CaCl2","pH_H2O","PointEC","N","P","K","OC","Moisture","Clay","Silt","Sand","CaCO3","Coarse"];

    if isempty(mesh)
        prior.hasMeshMapping = false;
        prior.cellIndex = [];
        prior.cellCenter = [];
        prior.xyDistance = [];
        return;
    end

    if ~isfield(mesh, 'elementCenters') || isempty(mesh.elementCenters)
        error('Mesh must contain elementCenters for prior-sample mapping.');
    end

    cellCenters = mesh.elementCenters;
    queryPts = [tbl.x_m, tbl.y_m, zeros(n, 1)];
    cellIndex = zeros(n, 1);
    cellCenter = zeros(n, 3);
    xyDistance = zeros(n, 1);

    for i = 1:n
        dxy = hypot(cellCenters(:, 1) - queryPts(i, 1), cellCenters(:, 2) - queryPts(i, 2));
        dz = abs(cellCenters(:, 3) - queryPts(i, 3));
        % Prefer shallow cells when multiple cells are similar in x-y.
        score = dxy + 0.15 * dz;
        [~, idx] = min(score);
        cellIndex(i) = idx;
        cellCenter(i, :) = cellCenters(idx, :);
        xyDistance(i) = dxy(idx);
    end

    prior.hasMeshMapping = true;
    prior.cellIndex = cellIndex;
    prior.cellCenter = cellCenter;
    prior.xyDistance = xyDistance;

    mapTbl = table( ...
        tbl.sample_id, tbl.x_m, tbl.y_m, cellIndex, ...
        cellCenter(:, 1), cellCenter(:, 2), cellCenter(:, 3), xyDistance, ...
        'VariableNames', ["sample_id","x_m","y_m","cell_index","cell_x","cell_y","cell_z","xy_distance"]);
    prior.mappingTable = mapTbl;
end
