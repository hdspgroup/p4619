function [deltaM, info] = solve_property_update(J, rData, Wd, R, rhsReg)
%SOLVE_PROPERTY_UPDATE Solve the augmented least-squares update.

    A = [Wd * J; R];
    b = [Wd * rData; rhsReg];

    H = A' * A;
    g = A' * b;
    diagH = full(diag(H));
    positiveDiag = diagH(isfinite(diagH) & diagH > 0);
    if isempty(positiveDiag)
        damping = 1e-8;
    else
        damping = max(1e-10 * median(positiveDiag), 1e-12);
    end
    deltaM = (H + damping * speye(size(H, 1))) \ g;

    info = struct();
    info.systemRows = size(A, 1);
    info.systemCols = size(A, 2);
    info.damping = damping;
    info.residualNorm = norm(A * deltaM - b);
end
