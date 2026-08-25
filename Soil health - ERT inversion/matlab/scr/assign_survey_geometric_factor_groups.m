function survey = assign_survey_geometric_factor_groups(survey, nGroups)
%ASSIGN_SURVEY_GEOMETRIC_FACTOR_GROUPS Group ABMN measurements by |K| quantiles.
%
% Groups are assigned using quantiles of abs(geometryFactorK), but ties are
% never split across groups. Group 1 corresponds to the largest |K| values.

    arguments
        survey struct
        nGroups (1,1) double {mustBeInteger, mustBePositive} = 3
    end

    if ~isfield(survey, 'geometryFactorK') || isempty(survey.geometryFactorK)
        error('Survey must contain geometryFactorK to assign groups.');
    end

    kAbs = abs(survey.geometryFactorK(:));
    if numel(kAbs) ~= survey.nMeasurements
        error('geometryFactorK length must match survey.nMeasurements.');
    end

    [uniqueKDesc, ~, mapToUnique] = unique(kAbs, 'sorted');
    uniqueKDesc = flipud(uniqueKDesc(:));
    nUnique = numel(uniqueKDesc);

    countsDesc = zeros(nUnique, 1);
    for i = 1:nUnique
        countsDesc(i) = nnz(kAbs == uniqueKDesc(i));
    end

    totalCount = numel(kAbs);
    targetCounts = (1:(nGroups - 1))' * (totalCount / nGroups);

    uniqueGroup = ones(nUnique, 1);
    groupIdx = 1;
    cumulative = 0;
    for i = 1:nUnique
        uniqueGroup(i) = groupIdx;
        cumulative = cumulative + countsDesc(i);
        remainingUnique = nUnique - i;
        remainingGroups = nGroups - groupIdx;
        if remainingGroups <= 0
            continue;
        end
        if cumulative >= targetCounts(groupIdx) && remainingUnique >= remainingGroups
            groupIdx = groupIdx + 1;
        end
    end

    measurementGroup = zeros(totalCount, 1);
    for i = 1:nUnique
        measurementGroup(kAbs == uniqueKDesc(i)) = uniqueGroup(i);
    end

    groupLabels = strings(nGroups, 1);
    groupRanges = zeros(nGroups, 2);
    for g = 1:nGroups
        mask = measurementGroup == g;
        if any(mask)
            kvals = kAbs(mask);
            groupRanges(g, :) = [min(kvals), max(kvals)];
            groupLabels(g) = sprintf('K-group %d | |K| in [%.3g, %.3g]', g, groupRanges(g, 1), groupRanges(g, 2));
        else
            groupRanges(g, :) = [NaN, NaN];
            groupLabels(g) = sprintf('K-group %d | empty', g);
        end
    end

    survey.grouping = struct();
    survey.grouping.type = "geometric_factor_quantile";
    survey.grouping.nGroups = nGroups;
    survey.grouping.groupIndex = measurementGroup;
    survey.grouping.groupLabels = groupLabels;
    survey.grouping.groupRanges = groupRanges;
    survey.grouping.kAbs = kAbs;
end
