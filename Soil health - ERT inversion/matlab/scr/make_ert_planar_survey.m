function survey = make_ert_planar_survey(options)
%MAKE_ERT_PLANAR_SURVEY Build a practical ABMN survey for arbitrary 2D electrode locations.
%
% The resulting survey is compatible with solve_ert_forward via make_ert_survey.

    arguments
        options.Electrodes (:,3) double
        options.Current (1,1) double = 1
        options.MinCurrentSpacing (1,1) double = 2.5
        options.MinPotentialSpacing (1,1) double = 2.5
        options.MinCentroidSeparation (1,1) double = 2.5
        options.MaxCurrentSpacing (1,1) double = inf
        options.MaxPotentialSpacing (1,1) double = inf
        options.MaxGeometricFactor (1,1) double = inf
        options.MinSensitivity (1,1) double = 0
        options.ShowProgress (1,1) logical = true
        options.ProgressMode (1,1) string = "auto" % auto | waitbar | text | off
        options.ProgressLabel (1,1) string = "Building practical ABMN combinations"
        options.TargetDepth (1,1) double = NaN
        options.DepthToCurrentSpacingFactor (1,1) double = 3.0
        options.DepthToPotentialSpacingFactor (1,1) double = 2.0
        options.DepthToCentroidSpacingFactor (1,1) double = 3.0
    end

    electrodes = options.Electrodes;
    nElectrodes = size(electrodes, 1);
    quads = zeros(0, 4);

    geometryFactorK = zeros(0, 1);
    geometryTerm = zeros(0, 1);
    distanceAB = zeros(0, 1);
    distanceMN = zeros(0, 1);
    centroidSeparation = zeros(0, 1);
    midpointAB = zeros(0, 2);
    midpointMN = zeros(0, 2);

    totalABPairs = nchoosek(nElectrodes, 2);
    processedABPairs = 0;
    totalCandidateQuads = 0;
    acceptedQuads = 0;
    lastPrintedPercent = -1;
    progressHandle = [];
    progressMode = localResolveProgressMode(options.ShowProgress, options.ProgressMode);
    useWaitbar = progressMode == "waitbar";

    maxCurrentSpacing = options.MaxCurrentSpacing;
    maxPotentialSpacing = options.MaxPotentialSpacing;
    minCentroidSeparation = options.MinCentroidSeparation;
    maxCentroidSpacing = inf;
    if isfinite(options.TargetDepth) && options.TargetDepth > 0
        maxCurrentSpacing = min(maxCurrentSpacing, options.DepthToCurrentSpacingFactor * options.TargetDepth);
        maxPotentialSpacing = min(maxPotentialSpacing, options.DepthToPotentialSpacingFactor * options.TargetDepth);
        maxCentroidSpacing = options.DepthToCentroidSpacingFactor * options.TargetDepth;
    end

    if options.ShowProgress
        if useWaitbar
            progressHandle = waitbar(0, char(options.ProgressLabel), ...
                'Name', 'ERT Survey Builder');
        else
            fprintf('%s\n', options.ProgressLabel);
        end
    end

    try
        for A = 1:(nElectrodes - 3)
            for B = (A + 1):(nElectrodes - 2)
                processedABPairs = processedABPairs + 1;
                updateProgress();

                dAB = norm(electrodes(A, 1:2) - electrodes(B, 1:2));
                if dAB < options.MinCurrentSpacing || dAB > maxCurrentSpacing
                    continue;
                end

                for M = 1:(nElectrodes - 1)
                    if any(M == [A, B])
                        continue;
                    end
                    for N = (M + 1):nElectrodes
                        if any(N == [A, B])
                            continue;
                        end

                        totalCandidateQuads = totalCandidateQuads + 1;

                        dMN = norm(electrodes(M, 1:2) - electrodes(N, 1:2));
                        if dMN < options.MinPotentialSpacing || dMN > maxPotentialSpacing
                            continue;
                        end

                        midAB = 0.5 * (electrodes(A, 1:2) + electrodes(B, 1:2));
                        midMN = 0.5 * (electrodes(M, 1:2) + electrodes(N, 1:2));
                        dCent = norm(midAB - midMN);
                        if dCent < minCentroidSeparation || dCent > maxCentroidSpacing
                            continue;
                        end

                        [K, g] = compute_ert_geometric_factor(electrodes, [A, B, M, N]);
                        if ~isfinite(K) || ~isfinite(g)
                            continue;
                        end
                        if abs(K) > options.MaxGeometricFactor
                            continue;
                        end
                        if abs(g) < options.MinSensitivity
                            continue;
                        end

                        quads(end+1, :) = [A, B, M, N]; %#ok<AGROW>
                        geometryFactorK(end+1, 1) = K; %#ok<AGROW>
                        geometryTerm(end+1, 1) = g; %#ok<AGROW>
                        distanceAB(end+1, 1) = dAB; %#ok<AGROW>
                        distanceMN(end+1, 1) = dMN; %#ok<AGROW>
                        centroidSeparation(end+1, 1) = dCent; %#ok<AGROW>
                        midpointAB(end+1, :) = midAB; %#ok<AGROW>
                        midpointMN(end+1, :) = midMN; %#ok<AGROW>
                        acceptedQuads = acceptedQuads + 1;
                    end
                end
            end
        end
    catch ME
        closeProgress();
        rethrow(ME);
    end

    if isempty(quads)
        error('No practical ABMN combinations passed the filters.');
    end

    survey = make_ert_survey('ertline', ...
        'Electrodes', electrodes, ...
        'Quads', quads, ...
        'Current', options.Current);

    survey.arrayLabel = repmat("all-practical-planar", survey.nMeasurements, 1);
    survey.arrayTypes = "all-practical-planar";
    survey.geometryFactorK = geometryFactorK;
    survey.geometryTerm = geometryTerm;
    survey.distanceAB = distanceAB;
    survey.distanceMN = distanceMN;
    survey.centroidSeparation = centroidSeparation;
    survey.midpointAB = midpointAB;
    survey.midpointMN = midpointMN;
    survey.apparentResistivityScale = geometryFactorK / options.Current;
    survey.electrodeUsage = localElectrodeUsage(quads, nElectrodes);
    survey.pairUsage = localPairUsage(quads, nElectrodes);
    survey.summary = struct( ...
        'totalABPairs', totalABPairs, ...
        'processedABPairs', processedABPairs, ...
        'totalCandidateQuads', totalCandidateQuads, ...
        'acceptedQuads', acceptedQuads, ...
        'acceptanceRate', acceptedQuads / max(totalCandidateQuads, 1), ...
        'unusedElectrodes', sum(survey.electrodeUsage.total == 0), ...
        'targetDepth', options.TargetDepth, ...
        'maxCurrentSpacingUsed', maxCurrentSpacing, ...
        'maxPotentialSpacingUsed', maxPotentialSpacing, ...
        'maxCentroidSpacingUsed', maxCentroidSpacing);

    closeProgress();
    printSummary();

    function updateProgress()
        if ~options.ShowProgress
            return;
        end

        fractionDone = processedABPairs / max(totalABPairs, 1);
        percentDone = floor(100 * fractionDone);

        if useWaitbar
            if ~isempty(progressHandle) && isgraphics(progressHandle)
                waitbar(fractionDone, progressHandle, sprintf('%s (%d/%d A-B pairs)', ...
                    options.ProgressLabel, processedABPairs, totalABPairs));
            end
            return;
        end

        if percentDone > lastPrintedPercent
            lastPrintedPercent = percentDone;
            fprintf('  %3d%% complete (%d/%d A-B pairs)\n', ...
                percentDone, processedABPairs, totalABPairs);
        end
    end

    function closeProgress()
        if ~options.ShowProgress
            return;
        end

        if useWaitbar
            if ~isempty(progressHandle) && isgraphics(progressHandle)
                close(progressHandle);
            end
        else
            fprintf('  100%% complete (%d/%d A-B pairs)\n', totalABPairs, totalABPairs);
        end
    end

    function printSummary()
        acceptanceRatePct = 100 * acceptedQuads / max(totalCandidateQuads, 1);
        fprintf('\nSurvey summary\n');
        fprintf('  A-B pairs checked: %d / %d\n', processedABPairs, totalABPairs);
        fprintf('  Candidate ABMN combinations: %d\n', totalCandidateQuads);
        fprintf('  Accepted ABMN combinations:  %d\n', acceptedQuads);
        fprintf('  Acceptance rate: %.2f%%\n', acceptanceRatePct);
        fprintf('  Unused electrodes: %d\n', sum(survey.electrodeUsage.total == 0));
        if isfinite(options.TargetDepth) && options.TargetDepth > 0
            fprintf('  Target depth: %.2f m\n', options.TargetDepth);
            fprintf('  Max current spacing used:   %.2f m\n', maxCurrentSpacing);
            fprintf('  Max potential spacing used: %.2f m\n', maxPotentialSpacing);
            fprintf('  Max centroid spacing used:  %.2f m\n', maxCentroidSpacing);
        end
        fprintf('\n');
    end
end

function usage = localElectrodeUsage(quads, nElectrodes)
    usage = struct();
    usage.A = accumarray(quads(:, 1), 1, [nElectrodes, 1], @sum, 0);
    usage.B = accumarray(quads(:, 2), 1, [nElectrodes, 1], @sum, 0);
    usage.M = accumarray(quads(:, 3), 1, [nElectrodes, 1], @sum, 0);
    usage.N = accumarray(quads(:, 4), 1, [nElectrodes, 1], @sum, 0);
    usage.current = usage.A + usage.B;
    usage.potential = usage.M + usage.N;
    usage.total = usage.current + usage.potential;
end

function pairUsage = localPairUsage(quads, nElectrodes)
    pairUsage = struct();
    pairUsage.current = zeros(nElectrodes, nElectrodes);
    pairUsage.potential = zeros(nElectrodes, nElectrodes);

    for i = 1:size(quads, 1)
        A = quads(i, 1); B = quads(i, 2);
        M = quads(i, 3); N = quads(i, 4);
        pairUsage.current(A, B) = pairUsage.current(A, B) + 1;
        pairUsage.current(B, A) = pairUsage.current(B, A) + 1;
        pairUsage.potential(M, N) = pairUsage.potential(M, N) + 1;
        pairUsage.potential(N, M) = pairUsage.potential(N, M) + 1;
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
