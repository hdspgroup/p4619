function info = evaluate_ec_background_anomaly_objective(dObs, dPred, Wd, mesh, param, uBg, uAnom, uBgRef, uAnomRef, weights)
%EVALUATE_EC_BACKGROUND_ANOMALY_OBJECTIVE Objective for background+anomaly parameterization.

    [Dx, Dy, Dz, Dxx, Dyy, Dzz] = make_property_difference_operators(mesh);
    [DzBg, DzzBg] = localBackgroundDifferenceOperators(param.zCenters);
    bgProfile = param.basisZ * uBg;

    dataResidual = Wd * (dObs - dPred);
    phiData = sum(dataResidual .^ 2);

    phiRegBg = weights.betaBg * sum((uBg - uBgRef) .^ 2) + ...
        weights.z1Bg * sum((DzBg * bgProfile) .^ 2) + ...
        weights.z2Bg * sum((DzzBg * bgProfile) .^ 2);

    phiRegAnom = weights.zeroMeanAnom * sum((param.Amean * uAnom) .^ 2) + ...
        weights.betaAnom * sum((uAnom - uAnomRef) .^ 2) + ...
        weights.axAnom * sum((Dx * uAnom) .^ 2) + ...
        weights.ayAnom * sum((Dy * uAnom) .^ 2) + ...
        weights.azAnom * sum((Dz * uAnom) .^ 2) + ...
        weights.cxAnom * sum((Dxx * uAnom) .^ 2) + ...
        weights.cyAnom * sum((Dyy * uAnom) .^ 2) + ...
        weights.czAnom * sum((Dzz * uAnom) .^ 2);

    info = struct();
    info.total = phiData + phiRegBg + phiRegAnom;
    info.phiData = phiData;
    info.phiReg = phiRegBg + phiRegAnom;
    info.phiRegBg = phiRegBg;
    info.phiRegAnom = phiRegAnom;
    info.weightedDataNorm = norm(dataResidual);
    info.normalizedDataMisfit = phiData / max(1, numel(dObs));
end

function [DzBg, DzzBg] = localBackgroundDifferenceOperators(zCenters)
    nZ = numel(zCenters);
    if nZ <= 1
        DzBg = sparse(0, nZ);
        DzzBg = sparse(0, nZ);
        return;
    end

    h = max(diff(zCenters(:)), eps);
    nRows1 = nZ - 1;
    rr = repelem((1:nRows1)', 2, 1);
    cc = [(1:nRows1)'; (2:nZ)'];
    vv = [-1 ./ h; 1 ./ h];
    DzBg = sparse(rr, cc, vv, nRows1, nZ);

    if nZ <= 2
        DzzBg = sparse(0, nZ);
    else
        nRows2 = nZ - 2;
        hLeft = h(1:end-1);
        hRight = h(2:end);
        denom = hLeft + hRight;
        rr = repelem((1:nRows2)', 3, 1);
        cc = [(1:nRows2)'; (2:nZ-1)'; (3:nZ)'];
        vv = [ ...
            2 ./ (hLeft .* denom); ...
            -2 .* (1 ./ hLeft + 1 ./ hRight) ./ denom; ...
            2 ./ (hRight .* denom)];
        DzzBg = sparse(rr, cc, vv, nRows2, nZ);
    end
end
