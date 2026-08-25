function info = evaluate_property_inversion_objective(m, dObs, dPred, Wd, config, mesh, mRef, weights, prior, sigmaCurrent)
%EVALUATE_PROPERTY_INVERSION_OBJECTIVE Evaluate the full inversion objective.

    if nargin < 9
        prior = [];
    end
    if nargin < 10
        sigmaCurrent = [];
    end

    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    fields = unpack_soil_fields(m, config);
    if isfield(weights, 'depthBasis')
        basisSettings = weights.depthBasis;
    else
        basisSettings = struct('nGaussians', 3, 'includeConstant', true);
    end
    Bdepth = make_depth_basis_matrix(mesh, basisSettings.nGaussians, basisSettings.includeConstant);
    Qdepth = orth(Bdepth);
    Pperp = sparse(speye(config.nCells) - Qdepth * Qdepth');

    dataResidual = Wd * (dObs - dPred);
    phiData = sum(dataResidual .^ 2);

    phiReg = 0;
    phiDepthBasis = 0;
    for i = 1:config.nProps
        name = config.names(i);
        idx = config.index.(char(name));
        u = m(idx);
        uRef = mRef(idx);
        duRef = u - uRef;
        w = weights.(char(name));
        phiReg = phiReg + w.beta * sum(duRef .^ 2);
        phiReg = phiReg + w.ax * sum((Dx * duRef) .^ 2);
        phiReg = phiReg + w.ay * sum((Dy * duRef) .^ 2);
        phiReg = phiReg + w.az * sum((Dz * duRef) .^ 2);
        phiReg = phiReg + w.cx * sum((Dxx * duRef) .^ 2);
        phiReg = phiReg + w.cy * sum((Dyy * duRef) .^ 2);
        phiReg = phiReg + w.cz * sum((Dzz * duRef) .^ 2);
        if isfield(w, 'depthBasisGamma') && w.depthBasisGamma > 0
            phiDepthBasis = phiDepthBasis + w.depthBasisGamma * sum((Pperp * duRef) .^ 2);
        end
    end

    textureResidual = fields.Clay + fields.Sand + fields.Silt - 100;
    phiTexture = weights.textureSumGamma * sum(textureResidual .^ 2);

    phiPrior = 0;
    phiPriorEC = 0;
    priorResidualSummary = struct();
    priorResidualSummary.usedProperties = strings(0, 1);
    priorResidualSummary.meanAbsResidual = struct();
    if ~isempty(prior) && isstruct(prior) && isfield(prior, 'hasMeshMapping') && prior.hasMeshMapping && ...
            isfield(weights, 'prior') && ~isempty(weights.prior)
        priorProps = string(fieldnames(weights.prior));
        used = strings(0, 1);
        for i = 1:numel(priorProps)
            name = priorProps(i);
            if ~ismember(name, config.names) || ~isfield(prior, char(name))
                continue;
            end
            gamma = weights.prior.(char(name));
            if ~(isnumeric(gamma) && isscalar(gamma) && gamma > 0)
                continue;
            end
            residual = fields.(char(name))(prior.cellIndex) - prior.(char(name));
            phiPrior = phiPrior + gamma * sum(residual .^ 2);
            used(end + 1, 1) = name; %#ok<AGROW>
            priorResidualSummary.meanAbsResidual.(char(name)) = mean(abs(residual));
        end
        priorResidualSummary.usedProperties = used;
    end

    if ~isempty(sigmaCurrent) && ~isempty(prior) && isstruct(prior) && isfield(prior, 'hasMeshMapping') && prior.hasMeshMapping && ...
            isfield(prior, 'PointEC') && isfield(weights, 'prior') && isfield(weights.prior, 'PointEC') && ...
            weights.prior.PointEC > 0
        ecResidual = sigmaCurrent(prior.cellIndex) - prior.PointEC;
        phiPriorEC = weights.prior.PointEC * sum(ecResidual .^ 2);
        priorResidualSummary.meanAbsResidual.PointEC = mean(abs(ecResidual));
    end

    info = struct();
    info.total = phiData + phiReg + phiDepthBasis + phiTexture + phiPrior + phiPriorEC;
    info.phiData = phiData;
    info.phiReg = phiReg;
    info.phiDepthBasis = phiDepthBasis;
    info.phiTexture = phiTexture;
    info.phiPrior = phiPrior;
    info.phiPriorEC = phiPriorEC;
    info.weightedDataNorm = norm(dataResidual);
    info.normalizedDataMisfit = phiData / max(1, numel(dObs));
    info.textureSumResidualMean = mean(abs(textureResidual));
    info.priorResidualSummary = priorResidualSummary;
end
