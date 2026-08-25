function A = assemble_ert_stiffness(mesh, sigmaCells, options)
%ASSEMBLE_ERT_STIFFNESS Assemble the finite-volume diffusion operator for ERT.
%
% The operator acts on cell-centered potentials and enforces local flux
% balance using two-point transmissibilities on an orthogonal rectilinear
% mesh. The function keeps the legacy name so the rest of the inversion
% stack can move over without a full vocabulary reset.

    arguments
        mesh struct
        sigmaCells (:,1) double
        options.UseParallel (1,1) logical = false %#ok<INUSA>
        options.ElementType = [] %#ok<INUSA>
        options.ShowProgress (1,1) logical = false
        options.ProgressMode (1,1) string = "auto"
    end

    if ~isfield(mesh, 'x') || ~isfield(mesh, 'y') || ~isfield(mesh, 'z')
        error('FVM assembly requires a rectilinear mesh with x, y, and z edge arrays.');
    end

    nx = numel(mesh.x) - 1;
    ny = numel(mesh.y) - 1;
    nz = numel(mesh.z) - 1;
    nCells = nx * ny * nz;

    if numel(sigmaCells) ~= nCells
        error('sigmaCells must contain one conductivity per cell (%d expected, %d received).', ...
            nCells, numel(sigmaCells));
    end

    progressHandle = [];
    progressMode = localResolveProgressMode(options.ShowProgress, options.ProgressMode);
    if options.ShowProgress
        if progressMode == "waitbar"
            progressHandle = waitbar(0, 'Assembling FVM ERT operator...');
        else
            fprintf('Assembling FVM ERT operator...\n');
        end
    end

    sigma3 = reshape(sigmaCells, [nx, ny, nz]);
    dx = diff(mesh.x(:)');
    dy = diff(mesh.y(:)');
    dz = diff(mesh.z(:)');

    idx3 = reshape(1:nCells, [nx, ny, nz]);
    diagVals = zeros(nCells, 1);
    I = [];
    J = [];
    S = [];

    updateProgress(0.1);

    % x-direction neighbors
    if nx > 1
        leftIdx = idx3(1:end-1, :, :);
        rightIdx = idx3(2:end, :, :);
        sigLeft = sigma3(1:end-1, :, :);
        sigRight = sigma3(2:end, :, :);
        areaX = reshape(dy, [1, ny, 1]) .* reshape(dz, [1, 1, nz]);
        distX = 0.5 * reshape(dx(1:end-1), [nx - 1, 1, 1]) + 0.5 * reshape(dx(2:end), [nx - 1, 1, 1]);
        Tx = localHarmonicMean(sigLeft, sigRight) .* areaX ./ distX;
        [diagVals, I, J, S] = localAccumulateDirection(diagVals, I, J, S, leftIdx, rightIdx, Tx);
    end
    updateProgress(0.4);

    % y-direction neighbors
    if ny > 1
        backIdx = idx3(:, 1:end-1, :);
        frontIdx = idx3(:, 2:end, :);
        sigBack = sigma3(:, 1:end-1, :);
        sigFront = sigma3(:, 2:end, :);
        areaY = reshape(dx, [nx, 1, 1]) .* reshape(dz, [1, 1, nz]);
        distY = 0.5 * reshape(dy(1:end-1), [1, ny - 1, 1]) + 0.5 * reshape(dy(2:end), [1, ny - 1, 1]);
        Ty = localHarmonicMean(sigBack, sigFront) .* areaY ./ distY;
        [diagVals, I, J, S] = localAccumulateDirection(diagVals, I, J, S, backIdx, frontIdx, Ty);
    end
    updateProgress(0.7);

    % z-direction neighbors
    if nz > 1
        topIdx = idx3(:, :, 1:end-1);
        bottomIdx = idx3(:, :, 2:end);
        sigTop = sigma3(:, :, 1:end-1);
        sigBottom = sigma3(:, :, 2:end);
        areaZ = reshape(dx, [nx, 1, 1]) .* reshape(dy, [1, ny, 1]);
        distZ = 0.5 * reshape(dz(1:end-1), [1, 1, nz - 1]) + 0.5 * reshape(dz(2:end), [1, 1, nz - 1]);
        Tz = localHarmonicMean(sigTop, sigBottom) .* areaZ ./ distZ;
        [diagVals, I, J, S] = localAccumulateDirection(diagVals, I, J, S, topIdx, bottomIdx, Tz);
    end
    updateProgress(0.95);

    A = sparse([I; (1:nCells)'], [J; (1:nCells)'], [S; diagVals], nCells, nCells);
    A = 0.5 * (A + A');
    updateProgress(1.0);
    closeProgress();

    function updateProgress(frac)
        if ~options.ShowProgress
            return;
        end
        if ~isempty(progressHandle) && isgraphics(progressHandle)
            waitbar(frac, progressHandle);
        elseif frac >= 1
            fprintf('  100%% complete\n');
        elseif frac >= 0.95
            fprintf('  95%% complete\n');
        elseif frac >= 0.70
            fprintf('  70%% complete\n');
        elseif frac >= 0.40
            fprintf('  40%% complete\n');
        elseif frac >= 0.10
            fprintf('  10%% complete\n');
        end
    end

    function closeProgress()
        if ~isempty(progressHandle) && isgraphics(progressHandle)
            close(progressHandle);
        end
    end
end

function [diagVals, I, J, S] = localAccumulateDirection(diagVals, I, J, S, idxA, idxB, T)
    ia = idxA(:);
    ib = idxB(:);
    tv = T(:);

    diagVals = diagVals + accumarray(ia, tv, size(diagVals), @sum, 0) ...
                        + accumarray(ib, tv, size(diagVals), @sum, 0);

    I = [I; ia; ib]; %#ok<AGROW>
    J = [J; ib; ia]; %#ok<AGROW>
    S = [S; -tv; -tv]; %#ok<AGROW>
end

function hm = localHarmonicMean(a, b)
    hm = 2 .* a .* b ./ max(a + b, realmin);
end

function progressMode = localResolveProgressMode(showProgress, requestedMode)
    if ~showProgress || strcmpi(requestedMode, "off")
        progressMode = "off";
        return;
    end

    requestedMode = lower(string(requestedMode));
    switch requestedMode
        case "auto"
            if usejava('desktop')
                progressMode = "waitbar";
            else
                progressMode = "text";
            end
        case {"waitbar", "text"}
            progressMode = requestedMode;
        otherwise
            error('Unsupported ProgressMode: %s', requestedMode);
    end
end
