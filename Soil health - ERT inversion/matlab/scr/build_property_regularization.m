function [R, rhsReg, info] = build_property_regularization(config, mesh, mCurrent, mRef, weights, prior)
%BUILD_PROPERTY_REGULARIZATION Assemble per-property roughness blocks in
% inversion-variable space, with texture-sum regularization linearized in
% physical space around the current model.

    if nargin < 6
        prior = [];
    end

    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    nCells = config.nCells;
    if isfield(weights, 'depthBasis')
        basisSettings = weights.depthBasis;
    else
        basisSettings = struct('nGaussians', 3, 'includeConstant', true);
    end
    Bdepth = make_depth_basis_matrix(mesh, basisSettings.nGaussians, basisSettings.includeConstant);
    Qdepth = orth(Bdepth);
    Pperp = sparse(speye(nCells) - Qdepth * Qdepth');

    Rblocks = cell(config.nProps, 1);
    rhsBlocks = cell(config.nProps, 1);

    for i = 1:config.nProps
        name = config.names(i);
        w = weights.(char(name));
        uCurrent = mCurrent(config.index.(char(name)));
        uRef = mRef(config.index.(char(name)));
        duRef = uCurrent - uRef;
        if isfield(w, 'depthBasisGamma') && w.depthBasisGamma > 0
            Rdepth = sqrt(w.depthBasisGamma) * Pperp;
            rhsDepth = -Rdepth * duRef;
        else
            Rdepth = sparse(0, nCells);
            rhsDepth = zeros(0, 1);
        end
        Ri = [ ...
            sqrt(w.beta) * speye(nCells); ...
            sqrt(w.ax) * Dx; ...
            sqrt(w.ay) * Dy; ...
            sqrt(w.az) * Dz; ...
            sqrt(w.cx) * Dxx; ...
            sqrt(w.cy) * Dyy; ...
            sqrt(w.cz) * Dzz; ...
            Rdepth];
        refBlock = -duRef;
        rhsi = [ ...
            sqrt(w.beta) * refBlock; ...
            -sqrt(w.ax) * (Dx * duRef); ...
            -sqrt(w.ay) * (Dy * duRef); ...
            -sqrt(w.az) * (Dz * duRef); ...
            -sqrt(w.cx) * (Dxx * duRef); ...
            -sqrt(w.cy) * (Dyy * duRef); ...
            -sqrt(w.cz) * (Dzz * duRef); ...
            rhsDepth];
        Rblocks{i} = Ri;
        rhsBlocks{i} = rhsi;
    end

    R = blkdiag(Rblocks{:});
    rhsReg = vertcat(rhsBlocks{:});

    [Rtex, rhsTex] = build_texture_sum_penalty(config, mCurrent, 100, weights.textureSumGamma);
    R = [R; Rtex];
    rhsReg = [rhsReg; rhsTex];

    [Rprior, rhsPrior, priorInfo] = build_property_prior_penalty(config, mCurrent, prior, weights);
    R = [R; Rprior];
    rhsReg = [rhsReg; rhsPrior];

    info = struct();
    info.nRowsBase = sum(cellfun(@(x) size(x, 1), Rblocks));
    info.nRowsTexture = size(Rtex, 1);
    info.nRowsPrior = size(Rprior, 1);
    info.nRowsTotal = size(R, 1);
    info.depthBasis = basisSettings;
    info.prior = priorInfo;
end
