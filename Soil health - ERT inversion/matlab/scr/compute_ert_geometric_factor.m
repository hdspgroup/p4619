function [K, g] = compute_ert_geometric_factor(electrodes, quad)
%COMPUTE_ERT_GEOMETRIC_FACTOR Geometric factor for a general ABMN quadrupole.
%
% For a half-space:
%   DeltaV = (rho * I / (2*pi)) * g
%   rho_a  = K * DeltaV / I
% with
%   K = 2*pi / g

    A = electrodes(quad(1), :);
    B = electrodes(quad(2), :);
    M = electrodes(quad(3), :);
    N = electrodes(quad(4), :);

    rAM = norm(A - M);
    rAN = norm(A - N);
    rBM = norm(B - M);
    rBN = norm(B - N);

    if any([rAM, rAN, rBM, rBN] <= 0)
        K = NaN;
        g = NaN;
        return;
    end

    g = 1 / rAM - 1 / rBM - 1 / rAN + 1 / rBN;
    if abs(g) < eps
        K = NaN;
    else
        K = 2 * pi / g;
    end
end
