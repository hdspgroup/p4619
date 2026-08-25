function info = evaluate_ec_inversion_objective(u, dObs, dPred, Wd, mesh, uRef, weights)
%EVALUATE_EC_INVERSION_OBJECTIVE Objective for EC-only inversion in log10 space.

    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);

    dataResidual = Wd * (dObs - dPred);
    phiData = sum(dataResidual .^ 2);
    phiReg = weights.beta * sum((u - uRef) .^ 2) + ...
        weights.ax * sum((Dx * u) .^ 2) + ...
        weights.ay * sum((Dy * u) .^ 2) + ...
        weights.az * sum((Dz * u) .^ 2) + ...
        weights.cx * sum((Dxx * u) .^ 2) + ...
        weights.cy * sum((Dyy * u) .^ 2) + ...
        weights.cz * sum((Dzz * u) .^ 2);

    info = struct();
    info.total = phiData + phiReg;
    info.phiData = phiData;
    info.phiReg = phiReg;
    info.weightedDataNorm = norm(dataResidual);
    info.normalizedDataMisfit = phiData / max(1, numel(dObs));
end

