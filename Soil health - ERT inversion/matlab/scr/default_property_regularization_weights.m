function weights = default_property_regularization_weights(csvPath)
%DEFAULT_PROPERTY_REGULARIZATION_WEIGHTS Per-property smoothing parameters.

    weights = struct();
    weights.depthBasis = struct('nGaussians', 3, 'includeConstant', true);
    weights.N = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.35, 'ax', 1.10, 'ay', 1.25, 'az', 1.60, ...
        'cx', 0.45, 'cy', 0.45, 'cz', 0.70, 'depthBasisGamma', 8.00);
    weights.P = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.30, 'ax', 0.85, 'ay', 1.00, 'az', 1.40, ...
        'cx', 0.35, 'cy', 0.35, 'cz', 0.58, 'depthBasisGamma', 7.00);
    weights.K = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.35, 'ax', 1.00, 'ay', 1.15, 'az', 1.55, ...
        'cx', 0.42, 'cy', 0.42, 'cz', 0.68, 'depthBasisGamma', 8.00);
    weights.OC = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.40, 'ax', 1.10, 'ay', 1.25, 'az', 1.70, ...
        'cx', 0.45, 'cy', 0.45, 'cz', 0.72, 'depthBasisGamma', 8.50);
    weights.Moisture = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.48, 'ax', 1.40, 'ay', 1.55, 'az', 2.00, ...
        'cx', 0.58, 'cy', 0.58, 'cz', 0.95, 'depthBasisGamma', 10.00);
    weights.Clay = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 1.20, 'ax', 3.80, 'ay', 4.10, 'az', 4.80, ...
        'cx', 1.80, 'cy', 1.80, 'cz', 2.40, 'depthBasisGamma', 22.00);
    weights.Sand = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.85, 'ax', 2.40, 'ay', 2.60, 'az', 3.20, ...
        'cx', 1.10, 'cy', 1.10, 'cz', 1.60, 'depthBasisGamma', 15.00);
    weights.Silt = struct('betaBg', 10, 'z1Bg', 10, 'z2Bg', 10, ...
        'beta', 0.85, 'ax', 2.40, 'ay', 2.60, 'az', 3.20, ...
        'cx', 1.10, 'cy', 1.10, 'cz', 1.60, 'depthBasisGamma', 15.00);
    weights.textureSumGamma = 2.50;
    weights.prior = struct( ...
        'PointEC', 2.00, ...
        'N', 0.80, ...
        'P', 0.40, ...
        'K', 0.40, ...
        'OC', 0.40, ...
        'Clay', 1.60, ...
        'Moisture', 0.90, ...
        'Silt', 1.00, ...
        'Sand', 1.00);

    if nargin >= 1 && ~isempty(csvPath) && isfile(csvPath)
        weights = localApplyWeightTable(weights, csvPath);
    end
end

function weights = localApplyWeightTable(weights, csvPath)
    raw = readcell(csvPath, 'FileType', 'text');
    raw = localTrimTrailingEmpty(raw);
    if isempty(raw)
        return;
    end

    % Left block: global_name/value pairs in columns A:B
    for i = 2:size(raw, 1)
        key = localCellString(raw{i, 1});
        if key == "" || any(strcmpi(key, ["background_name", "property"]))
            continue;
        end
        value = localCellNumber(raw{i, 2});
        if ~isfinite(value)
            continue;
        end
        switch lower(key)
            case "texturesumgamma"
                weights.textureSumGamma = value;
            case "priorpointec"
                weights.prior.PointEC = value;
            case "depthbasisngaussians"
                weights.depthBasis.nGaussians = value;
            case "depthbasisincludeconstant"
                weights.depthBasis.includeConstant = logical(value);
        end
    end

    % Right block: property table starting at column D
    headerRow = [];
    for i = 1:size(raw, 1)
        if size(raw, 2) >= 4 && strcmpi(localCellString(raw{i, 4}), "property")
            headerRow = i;
            break;
        end
    end
    if isempty(headerRow)
        return;
    end

    propHeaders = strings(1, size(raw, 2) - 3);
    for j = 4:size(raw, 2)
        propHeaders(j - 3) = localCellString(raw{headerRow, j});
    end

    propertyFieldNames = string(fieldnames(weights));
    for i = headerRow + 1:size(raw, 1)
        rowName = strtrim(localCellString(raw{i, 4}));
        if rowName == "" || strcmpi(rowName, "EC") || strcmpi(rowName, "EC_ANOMALY")
            continue;
        end
        fieldIdx = find(strcmpi(propertyFieldNames, rowName), 1);
        if isempty(fieldIdx)
            continue;
        end
        weightFieldName = char(propertyFieldNames(fieldIdx));

        fieldNames = ["betaBg","z1Bg","z2Bg","beta","ax","ay","az","cx","cy","cz","depthBasisGamma"];
        for j = 1:numel(fieldNames)
            fieldName = fieldNames(j);
            colIdx = localFindWeightHeader(propHeaders, fieldName);
            if isempty(colIdx)
                continue;
            end
            value = localCellNumber(raw{i, colIdx + 3});
            if isfinite(value)
                weights.(weightFieldName).(char(fieldName)) = value;
            end
        end

        colIdx = localFindWeightHeader(propHeaders, "priorGamma");
        priorFieldNames = string(fieldnames(weights.prior));
        priorIdx = find(strcmpi(priorFieldNames, rowName), 1);
        if ~isempty(colIdx) && ~isempty(priorIdx)
            value = localCellNumber(raw{i, colIdx + 3});
            if isfinite(value)
                weights.prior.(char(priorFieldNames(priorIdx))) = value;
            end
        end
    end
end

function colIdx = localFindWeightHeader(headers, fieldName)
    aliases = localWeightHeaderAliases(fieldName);
    colIdx = [];
    for i = 1:numel(aliases)
        colIdx = find(strcmpi(headers, aliases(i)), 1);
        if ~isempty(colIdx)
            return;
        end
    end
end

function aliases = localWeightHeaderAliases(fieldName)
    switch string(fieldName)
        case "z1Bg"
            aliases = ["z1Bg", "z1Bg_grad_per_m"];
        case "z2Bg"
            aliases = ["z2Bg", "z2Bg_curv_per_m2"];
        case "ax"
            aliases = ["ax", "ax_grad_per_m"];
        case "ay"
            aliases = ["ay", "ay_grad_per_m"];
        case "az"
            aliases = ["az", "az_grad_per_m"];
        case "cx"
            aliases = ["cx", "cx_curv_per_m2"];
        case "cy"
            aliases = ["cy", "cy_curv_per_m2"];
        case "cz"
            aliases = ["cz", "cz_curv_per_m2"];
        otherwise
            aliases = string(fieldName);
    end
end

function raw = localTrimTrailingEmpty(raw)
    if isempty(raw)
        return;
    end
    keepRows = false(size(raw, 1), 1);
    for i = 1:size(raw, 1)
        for j = 1:size(raw, 2)
            if localCellString(raw{i, j}) ~= "" || (isnumeric(raw{i, j}) && ~isempty(raw{i, j}) && isfinite(raw{i, j}))
                keepRows(i) = true;
            end
        end
    end
    raw = raw(keepRows, :);

    keepUntil = size(raw, 2);
    while keepUntil > 1
        colHasValue = false;
        for i = 1:size(raw, 1)
            v = raw{i, keepUntil};
            if localCellString(v) ~= "" || (isnumeric(v) && ~isempty(v) && isfinite(v))
                colHasValue = true;
                break;
            end
        end
        if colHasValue
            break;
        end
        keepUntil = keepUntil - 1;
    end
    raw = raw(:, 1:keepUntil);
end

function s = localCellString(v)
    if isempty(v)
        s = "";
    elseif iscell(v)
        if isempty(v)
            s = "";
        else
            s = localCellString(v{1});
        end
    elseif isstring(v)
        if isempty(v)
            s = "";
        else
            s = strtrim(v(1));
        end
    elseif ischar(v)
        s = strtrim(string(v));
    elseif isnumeric(v)
        s = "";
    elseif islogical(v)
        s = "";
    else
        try
            if all(ismissing(v))
                s = "";
            else
                sv = string(v);
                if isempty(sv)
                    s = "";
                else
                    s = strtrim(sv(1));
                end
            end
        catch
            s = "";
        end
    end
end

function x = localCellNumber(v)
    if isnumeric(v) && isscalar(v)
        x = v;
        return;
    end
    s = localCellString(v);
    if s == ""
        x = NaN;
    else
        x = str2double(s);
    end
end
