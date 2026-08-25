function B = make_depth_basis_matrix(mesh, nGaussians, includeConstant)
%MAKE_DEPTH_BASIS_MATRIX Build broad depth-only basis functions over cells.
%
% B is [nCells x nBasis]. Each column depends only on cell-center depth.

    if nargin < 2 || isempty(nGaussians)
        nGaussians = 3;
    end
    if nargin < 3
        includeConstant = true;
    end

    z = mesh.elementCenters(:, 3);
    zMin = min(z);
    zMax = max(z);

    cols = {};
    if includeConstant
        cols{end + 1} = ones(size(z)); %#ok<AGROW>
    end

    if nGaussians > 0
        centers = linspace(zMin, zMax, nGaussians + 2);
        centers = centers(2:end-1);
        if numel(centers) > 1
            spacing = mean(diff(centers));
        else
            spacing = max((zMax - zMin) / 2, eps);
        end
        sigma = max(0.85 * spacing, max((zMax - zMin) / 6, eps));

        for i = 1:numel(centers)
            cols{end + 1} = exp(-0.5 * ((z - centers(i)) ./ sigma) .^ 2); %#ok<AGROW>
        end
    end

    B = zeros(numel(z), numel(cols));
    for i = 1:numel(cols)
        B(:, i) = cols{i};
    end

    % Normalize columns for better conditioning.
    for i = 1:size(B, 2)
        ni = norm(B(:, i));
        if ni > 0
            B(:, i) = B(:, i) / ni;
        end
    end
end
