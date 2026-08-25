function [Rprior, rhsPrior, info] = build_property_prior_penalty(config, mCurrent, prior, weights)
%BUILD_PROPERTY_PRIOR_PENALTY Assemble soft prior rows from sample data.
%
% The prior file is interpreted in physical property units. This penalty is
% linearized around the current inversion model so it fits the existing
% least-squares update:
%
%   x(u + du) ~= x(u) + dx/du * du
%
% and softly enforces:
%
%   x(sampleCell) ~= x_prior

    if nargin < 4 || isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'hasMeshMapping') || ~prior.hasMeshMapping
        Rprior = sparse(0, config.totalSize);
        rhsPrior = zeros(0, 1);
        info = struct('nRows', 0, 'usedProperties', strings(0, 1), 'nUnusedSamples', 0);
        return;
    end

    if ~isfield(weights, 'prior') || isempty(weights.prior)
        Rprior = sparse(0, config.totalSize);
        rhsPrior = zeros(0, 1);
        info = struct('nRows', 0, 'usedProperties', strings(0, 1), 'nUnusedSamples', 0);
        return;
    end

    fields = unpack_soil_fields(mCurrent, config);
    propNames = string(fieldnames(weights.prior));

    rowCount = 0;
    usedProperties = strings(0, 1);
    for i = 1:numel(propNames)
        name = propNames(i);
        if ~ismember(name, config.names) || ~isfield(prior, char(name))
            continue;
        end
        gamma = weights.prior.(char(name));
        if ~(isnumeric(gamma) && isscalar(gamma) && gamma > 0)
            continue;
        end
        rowCount = rowCount + prior.nSamples;
        usedProperties(end + 1, 1) = name; %#ok<AGROW>
    end

    if rowCount == 0
        Rprior = sparse(0, config.totalSize);
        rhsPrior = zeros(0, 1);
        info = struct('nRows', 0, 'usedProperties', strings(0, 1), 'nUnusedSamples', prior.nSamples);
        return;
    end

    rr = zeros(rowCount, 1);
    cc = zeros(rowCount, 1);
    vv = zeros(rowCount, 1);
    rhsPrior = zeros(rowCount, 1);
    row = 0;

    for i = 1:numel(usedProperties)
        name = usedProperties(i);
        gamma = weights.prior.(char(name));
        idxAll = config.index.(char(name));
        fieldsCurrent = fields.(char(name));
        uCurrent = mCurrent(idxAll);
        dxdu = d_soil_property_d_inversion_var(name, uCurrent, config);
        targetVals = prior.(char(name));

        for j = 1:prior.nSamples
            cellIdx = prior.cellIndex(j);
            row = row + 1;
            rr(row) = row;
            cc(row) = idxAll(cellIdx);
            vv(row) = sqrt(gamma) * dxdu(cellIdx);
            rhsPrior(row) = sqrt(gamma) * (targetVals(j) - fieldsCurrent(cellIdx));
        end
    end

    Rprior = sparse(rr, cc, vv, rowCount, config.totalSize);

    info = struct();
    info.nRows = rowCount;
    info.usedProperties = usedProperties;
    info.nUnusedSamples = 0;
end
