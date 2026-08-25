function Ju = build_ert_jacobian_fd_u(forwardMesh, survey, uBase, fineFromCoarse, baseResult, options)
%BUILD_ERT_JACOBIAN_FD_U One-sided finite-difference Jacobian for coarse log10(EC).
%
% The unknowns live on the coarse inversion mesh in u = log10(EC). Each
% perturbation is mapped to the fine forward mesh before solving the
% electric potential field.

    arguments
        forwardMesh struct
        survey struct
        uBase (:,1) double
        fineFromCoarse double
        baseResult struct
        options.UseParallel (1,1) logical = false
        options.ShowProgress (1,1) logical = true
        options.ProgressMode (1,1) string = "auto" % auto | waitbar | text | off
        options.RelativePerturbation (1,1) double = 0.02
        options.AbsolutePerturbation (1,1) double = 1e-3
    end

    nParams = numel(uBase);
    nData = numel(baseResult.voltage);
    Ju = zeros(nData, nParams);

    progressCount = 0;
    progressLastPrint = 0;
    progressHandle = [];
    progressMode = localResolveProgressMode(options.ShowProgress, options.ProgressMode);
    if options.ShowProgress
        if progressMode == "waitbar"
            progressHandle = waitbar(0, 'Building one-sided Jacobian in log10(EC)...');
        else
            fprintf('Building one-sided Jacobian in log10(EC)...\n');
        end
    end

    if options.UseParallel
        pool = gcp('nocreate');
        if isempty(pool)
            parpool;
            pool = gcp('nocreate');
        end
        if options.ShowProgress
            fprintf('  J_u setup: %d parameters, %d data, parallel over columns with %d workers\n', ...
                nParams, nData, pool.NumWorkers);
        end
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) updateProgress());

        cols = cell(nParams, 1);
        parfor i = 1:nParams
            cols{i} = localForwardDifferenceColumn( ...
                i, forwardMesh, survey, uBase, fineFromCoarse, baseResult, ...
                options.RelativePerturbation, options.AbsolutePerturbation);
            send(q, 1);
        end
        for i = 1:nParams
            Ju(:, i) = cols{i};
        end
    else
        if options.ShowProgress
            fprintf('  J_u setup: %d parameters, %d data, serial over columns\n', nParams, nData);
        end
        for i = 1:nParams
            Ju(:, i) = localForwardDifferenceColumn( ...
                i, forwardMesh, survey, uBase, fineFromCoarse, baseResult, ...
                options.RelativePerturbation, options.AbsolutePerturbation);
            updateProgress();
        end
    end

    closeProgress();

    function updateProgress()
        if ~options.ShowProgress
            return;
        end
        progressCount = progressCount + 1;
        frac = progressCount / nParams;
        if ~isempty(progressHandle) && isgraphics(progressHandle)
            waitbar(frac, progressHandle);
        elseif frac - progressLastPrint >= 0.05 || frac == 1
            fprintf('  J_u %3.0f%%\n', 100 * frac);
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

function col = localForwardDifferenceColumn(iParam, forwardMesh, survey, uBase, fineFromCoarse, baseResult, relativePerturbation, absolutePerturbation)
%LOCALFORWARDDIFFERENCECOLUMN One Jacobian column in coarse log10(EC).

    du = max(absolutePerturbation, relativePerturbation * max(abs(uBase(iParam)), 1));

    uPert = uBase;
    uPert(iParam) = uPert(iParam) + du;
    sigmaPert = 10 .^ (fineFromCoarse * uPert);

    Kpert = assemble_ert_stiffness(forwardMesh, sigmaPert, ...
        'UseParallel', false, ...
        'ElementType', forwardMesh.elementType, ...
        'ShowProgress', false);
    resultPert = solve_ert_forward(forwardMesh, survey, Kpert, ...
        'ReferenceNode', baseResult.referenceNode);

    col = (resultPert.voltage - baseResult.voltage) / du;
end
