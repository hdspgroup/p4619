function sub = subset_ert_survey(survey, mask)
%SUBSET_ERT_SURVEY Return a measurement subset of an ERT survey struct.

    arguments
        survey struct
        mask (:,1) logical
    end

    if numel(mask) ~= survey.nMeasurements
        error('Mask length must match survey.nMeasurements.');
    end

    sub = survey;
    sub.quads = survey.quads(mask, :);
    sub.nMeasurements = size(sub.quads, 1);

    measurementFields = { ...
        'arrayLabel', 'geometryFactorK', 'geometryTerm', 'distanceAB', ...
        'distanceMN', 'centroidSeparation', 'midpointAB', 'midpointMN', ...
        'apparentResistivityScale'};
    for i = 1:numel(measurementFields)
        f = measurementFields{i};
        if isfield(survey, f) && ~isempty(survey.(f))
            val = survey.(f);
            if size(val, 1) == survey.nMeasurements
                sub.(f) = val(mask, :);
            elseif isvector(val) && numel(val) == survey.nMeasurements
                sub.(f) = val(mask);
            end
        end
    end

    if isfield(survey, 'grouping') && isfield(survey.grouping, 'groupIndex')
        sub.grouping = survey.grouping;
        sub.grouping.groupIndex = survey.grouping.groupIndex(mask);
        sub.grouping.kAbs = survey.grouping.kAbs(mask);
    end
end
