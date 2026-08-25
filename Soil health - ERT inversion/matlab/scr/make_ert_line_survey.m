function survey = make_ert_line_survey(options)
%MAKE_ERT_LINE_SURVEY Build a line survey and quadrupole list.
%
% Supported array types:
%   "wenner"
%   "dipole-dipole"
%   "schlumberger"
%   "all-practical"
%
% ArrayType can be a scalar string or a string array. When multiple array
% types are requested, the resulting quadrupoles are concatenated and
% duplicate rows are removed.

    arguments
        options.XElectrodes double = -25:5:25
        options.Y (1,1) double = 0
        options.Z (1,1) double = 0
        options.Current (1,1) double = 1
        options.ArrayType {mustBeTextOrTextArray} = "wenner"
        options.MinCurrentDipole (1,1) double = 1
        options.MinPotentialDipole (1,1) double = 1
        options.MinDipoleGap (1,1) double = 1
        options.MaxSpanElectrodes (1,1) double = inf
        options.MinGeometryValue (1,1) double = 0
    end

    x = options.XElectrodes(:);
    electrodes = [x, options.Y * ones(numel(x), 1), options.Z * ones(numel(x), 1)];
    nElectrodes = numel(x);

    arrayTypes = string(options.ArrayType);
    quads = zeros(0, 4);
    arrayLabel = strings(0, 1);

    for k = 1:numel(arrayTypes)
        thisType = lower(strtrim(arrayTypes(k)));
        switch thisType
            case "wenner"
                q = localWennerQuads(nElectrodes);
            case "dipole-dipole"
                q = localDipoleDipoleQuads(nElectrodes);
            case "schlumberger"
                q = localSchlumbergerQuads(nElectrodes);
            case "all-practical"
                q = localAllPracticalQuads(x, options);
            otherwise
                error('Unsupported array type: %s', arrayTypes(k));
        end

        quads = [quads; q]; %#ok<AGROW>
        arrayLabel = [arrayLabel; repmat(thisType, size(q, 1), 1)]; %#ok<AGROW>
    end

    if isempty(quads)
        error('The requested survey configuration produced zero quadrupoles.');
    end

    [quadsUnique, ia] = unique(quads, 'rows', 'stable');
    arrayLabel = arrayLabel(ia);

    survey = make_ert_survey('ertline', ...
        'Electrodes', electrodes, ...
        'Quads', quadsUnique, ...
        'Current', options.Current);
    survey.arrayLabel = arrayLabel;
    survey.arrayTypes = unique(arrayLabel, 'stable');
end

function quads = localWennerQuads(nElectrodes)
    quads = zeros(0, 4);
    for a = 1:floor((nElectrodes - 1) / 3)
        for s = 1:(nElectrodes - 3 * a)
            quads(end+1, :) = [s, s + 3 * a, s + a, s + 2 * a]; %#ok<AGROW>
        end
    end
end

function quads = localDipoleDipoleQuads(nElectrodes)
    quads = zeros(0, 4);
    maxDipole = floor((nElectrodes - 1) / 3);
    for a = 1:maxDipole
        maxN = floor((nElectrodes - 1) / a) - 2;
        for n = 1:maxN
            maxStart = nElectrodes - (n + 2) * a;
            for s = 1:maxStart
                quads(end+1, :) = [s, s + a, s + (n + 1) * a, s + (n + 2) * a]; %#ok<AGROW>
            end
        end
    end
end

function quads = localSchlumbergerQuads(nElectrodes)
    quads = zeros(0, 4);
    maxCurrentHalfSpan = floor((nElectrodes - 1) / 2);
    for c = 2:maxCurrentHalfSpan
        for p = 1:(c - 1)
            maxStart = nElectrodes - 2 * c;
            for s = 1:maxStart
                quads(end+1, :) = [s, s + 2 * c, s + (c - p), s + (c + p)]; %#ok<AGROW>
            end
        end
    end
end

function quads = localAllPracticalQuads(xElectrodes, options)
    nElectrodes = numel(xElectrodes);
    quads = zeros(0, 4);

    for A = 1:(nElectrodes - 3)
        for B = (A + 1):(nElectrodes - 2)
            currentDipole = B - A;
            if currentDipole < options.MinCurrentDipole
                continue;
            end

            for M = (B + 1):(nElectrodes - 1)
                gap = M - B;
                if gap < options.MinDipoleGap
                    continue;
                end

                for N = (M + 1):nElectrodes
                    potentialDipole = N - M;
                    if potentialDipole < options.MinPotentialDipole
                        continue;
                    end

                    span = N - A;
                    if span > options.MaxSpanElectrodes
                        continue;
                    end

                    geometryValue = abs(localGeometricSensitivityValue(xElectrodes, [A, B, M, N]));
                    if geometryValue <= options.MinGeometryValue
                        continue;
                    end

                    quads(end+1, :) = [A, B, M, N]; %#ok<AGROW>
                end
            end
        end
    end
end

function g = localGeometricSensitivityValue(xElectrodes, ids)
%LOCALGEOMETRICSENSITIVITYVALUE Line-electrode geometry term used as a weak-measurement filter.

    x = xElectrodes(ids);
    rAM = abs(x(1) - x(3));
    rAN = abs(x(1) - x(4));
    rBM = abs(x(2) - x(3));
    rBN = abs(x(2) - x(4));

    g = 1 ./ rAM - 1 ./ rAN - 1 ./ rBM + 1 ./ rBN;
end

function mustBeTextOrTextArray(value)
    if ~(isstring(value) || ischar(value) || iscellstr(value))
        error('ArrayType must be a string, char vector, or cell array of character vectors.');
    end
end
