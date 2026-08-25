function Jert = build_ert_jacobian_fd(mesh, survey, sigmaBase, baseResult, options)
%BUILD_ERT_JACOBIAN_FD Finite-difference ERT Jacobian with respect to sigma.

    arguments
        mesh struct
        survey struct
        sigmaBase (:,1) double
        baseResult struct
        options.UseParallel (1,1) logical = false
        options.ShowProgress (1,1) logical = true
        options.RelativePerturbation (1,1) double = 0.02
    end

    nCells = numel(sigmaBase);
    nData = numel(baseResult.voltage);
    Jert = zeros(nData, nCells);

    progressCount = 0;
    progressLastPrint = 0;
    progressHandle = [];
    if options.ShowProgress
        if usejava('desktop')
            progressHandle = waitbar(0, 'Building ERT Jacobian by finite differences...');
        else
            fprintf('Building ERT Jacobian by finite differences...\n');
        end
    end

    if options.UseParallel
        pool = gcp('nocreate');
        if isempty(pool)
            parpool;
        end
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) updateProgress());

        Jcols = cell(nCells, 1);
        parfor i = 1:nCells
            Jcols{i} = build_ert_jacobian_fd_column( ...
                i, mesh, survey, sigmaBase, baseResult, options.RelativePerturbation);
            send(q, 1);
        end
        for i = 1:nCells
            Jert(:, i) = Jcols{i};
        end
    else
        for i = 1:nCells
            Jert(:, i) = build_ert_jacobian_fd_column( ...
                i, mesh, survey, sigmaBase, baseResult, options.RelativePerturbation);
            updateProgress();
        end
    end

    closeProgress();

    function updateProgress()
        if ~options.ShowProgress
            return;
        end
        progressCount = progressCount + 1;
        frac = progressCount / nCells;
        if ~isempty(progressHandle) && isgraphics(progressHandle)
            waitbar(frac, progressHandle);
        elseif frac - progressLastPrint >= 0.05 || frac == 1
            fprintf('  J_ERT %3.0f%%\n', 100 * frac);
            progressLastPrint = frac;
        end
    end

    function closeProgress()
        if ~isempty(progressHandle) && isgraphics(progressHandle)
            close(progressHandle);
            progressHandle = [];
        end
    end
end

function col = build_ert_jacobian_fd_column(iCell, mesh, survey, sigmaBase, baseResult, relativePerturbation)
%BUILD_ERT_JACOBIAN_FD_COLUMN One finite-difference Jacobian column.

    sigmaPert = sigmaBase;
    ds = max(1e-6, relativePerturbation * max(sigmaBase(iCell), 1e-6));
    sigmaPert(iCell) = sigmaPert(iCell) + ds;

    Kpert = assemble_ert_stiffness(mesh, sigmaPert, ...
        'UseParallel', false, ...
        'ElementType', mesh.elementType, ...
        'ShowProgress', false);
    resultPert = solve_ert_forward(mesh, survey, Kpert, ...
        'ReferenceNode', baseResult.referenceNode);
    col = (resultPert.voltage - baseResult.voltage) / ds;
end
