function weights = default_ec_background_anomaly_weights(csvPath)
%DEFAULT_EC_BACKGROUND_ANOMALY_WEIGHTS Regularization for background+anomaly EC inversion.

    weights = struct();

    % Depth-dependent background mode u_bg(z).
    weights.backgroundMode = "gaussian3";
    weights.betaBg = 10^1.5;
    weights.z1Bg = 10^1.5;
    weights.z2Bg = 10^1.0;

    % Cellwise anomaly mode u_anom(x,z).
    weights.zeroMeanAnom = 10^3.0;
    weights.betaAnom = 10^3.0;
    weights.axAnom = 10^2.5;
    weights.ayAnom = 10^2.5;
    weights.azAnom = 10^2.0;
    weights.cxAnom = 10^2.5;
    weights.cyAnom = 10^2.5;
    weights.czAnom = 10^2.0;

    if nargin >= 1 && ~isempty(csvPath) && isfile(csvPath)
        weights = localApplyEcWeightTable(weights, csvPath);
    end
end

function weights = localApplyEcWeightTable(weights, csvPath)
    raw = readcell(csvPath, 'FileType', 'text');
    raw = localTrimTrailingEmpty(raw);
    if isempty(raw)
        return;
    end

    % Legacy left-block background settings in columns A:B.
    for i = 2:size(raw, 1)
        key = upper(strtrim(localCellString(raw{i, 1})));
        if key == "ECBACKGROUNDMODE"
            valueText = strtrim(localCellString(raw{i, 2}));
            if valueText ~= ""
                weights.backgroundMode = string(valueText);
            end
        elseif any(key == ["BETABG","Z1BG","Z2BG"])
            value = localCellNumber(raw{i, 2});
            if ~isfinite(value)
                continue;
            end
            switch key
                case "BETABG"
                    weights.betaBg = value;
                case "Z1BG"
                    weights.z1Bg = value;
                case "Z2BG"
                    weights.z2Bg = value;
            end
        end
    end

    % Unified right-block table starting at column D.
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

    for i = headerRow + 1:size(raw, 1)
        rowName = upper(strtrim(localCellString(raw{i, 4})));
        if rowName == "EC" || rowName == "EC_ANOMALY"
            weights = localApplyEcRow(weights, raw(i, :), propHeaders);
        end
    end
end

function weights = localApplyEcRow(weights, rowCells, headers)
    bgModeIdx = localFindWeightHeader(headers, "backgroundMode");
    if ~isempty(bgModeIdx)
        valueText = strtrim(localCellString(rowCells{1, bgModeIdx + 3}));
        if valueText ~= ""
            weights.backgroundMode = string(valueText);
        end
    end
    weights.betaBg = localAssignMapped(headers, rowCells, "betaBg", weights.betaBg);
    weights.z1Bg = localAssignMapped(headers, rowCells, "z1Bg", weights.z1Bg);
    weights.z2Bg = localAssignMapped(headers, rowCells, "z2Bg", weights.z2Bg);
    weights.zeroMeanAnom = localAssignMapped(headers, rowCells, "zeroMean", weights.zeroMeanAnom);
    weights.betaAnom = localAssignMapped(headers, rowCells, "beta", weights.betaAnom);
    weights.axAnom = localAssignMapped(headers, rowCells, "ax", weights.axAnom);
    weights.ayAnom = localAssignMapped(headers, rowCells, "ay", weights.ayAnom);
    weights.azAnom = localAssignMapped(headers, rowCells, "az", weights.azAnom);
    weights.cxAnom = localAssignMapped(headers, rowCells, "cx", weights.cxAnom);
    weights.cyAnom = localAssignMapped(headers, rowCells, "cy", weights.cyAnom);
    weights.czAnom = localAssignMapped(headers, rowCells, "cz", weights.czAnom);
end

function valueOut = localAssignMapped(headers, rowCells, headerName, valueIn)
    valueOut = valueIn;
    colIdx = localFindWeightHeader(headers, headerName);
    if isempty(colIdx)
        return;
    end
    value = localCellNumber(rowCells{1, colIdx + 3});
    if isfinite(value)
        valueOut = value;
    end
end

function colIdx = localFindWeightHeader(headers, headerName)
    aliases = localWeightHeaderAliases(headerName);
    colIdx = [];
    for i = 1:numel(aliases)
        colIdx = find(strcmpi(headers, aliases(i)), 1);
        if ~isempty(colIdx)
            return;
        end
    end
end

function aliases = localWeightHeaderAliases(headerName)
    switch string(headerName)
        case "backgroundMode"
            aliases = ["backgroundMode", "ecBackgroundMode", "bgMode"];
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
            aliases = string(headerName);
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
