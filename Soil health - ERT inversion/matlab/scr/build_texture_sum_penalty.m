function [Rtex, rhsTex] = build_texture_sum_penalty(config, mCurrent, targetTotal, gamma)
%BUILD_TEXTURE_SUM_PENALTY Softly enforce Clay + Sand + Silt = targetTotal.
%
% The inversion vector can be transformed, so this term is linearized in
% physical-space texture fields about the current model.

    nCells = config.nCells;
    clayCols = config.index.Clay;
    sandCols = config.index.Sand;
    siltCols = config.index.Silt;

    fields = unpack_soil_fields(mCurrent, config);
    uClay = mCurrent(clayCols);
    uSand = mCurrent(sandCols);
    uSilt = mCurrent(siltCols);
    dClay = d_soil_property_d_inversion_var("Clay", uClay, config);
    dSand = d_soil_property_d_inversion_var("Sand", uSand, config);
    dSilt = d_soil_property_d_inversion_var("Silt", uSilt, config);

    rr = repelem((1:nCells)', 3, 1);
    cc = [clayCols(:); sandCols(:); siltCols(:)];
    vv = sqrt(gamma) * [dClay(:); dSand(:); dSilt(:)];
    Rtex = sparse(rr, cc, vv, nCells, config.totalSize);

    residual = targetTotal - (fields.Clay + fields.Sand + fields.Silt);
    rhsTex = sqrt(gamma) * residual;
end
