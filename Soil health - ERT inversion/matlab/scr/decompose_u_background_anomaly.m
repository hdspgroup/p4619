function [uBg, uAnom] = decompose_u_background_anomaly(mesh, uCell, param)
%DECOMPOSE_U_BACKGROUND_ANOMALY Split a cellwise log10(EC) field into background + anomaly.

    if nargin < 3 || isempty(param)
        param = make_background_anomaly_parameterization(mesh);
    end

    n = mesh.gridSize - 1;
    u3 = reshape(uCell, n);
    uBgProfile = squeeze(mean(mean(u3, 1), 2));
    uBgProfile = uBgProfile(:);

    uBg = param.basisZ \ uBgProfile;
    B = param.B;
    uAnom = uCell - B * uBg;
end
