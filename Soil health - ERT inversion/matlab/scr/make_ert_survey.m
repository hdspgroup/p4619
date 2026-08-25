function survey = make_ert_survey(type, options)
%MAKE_ERT_SURVEY Create a generic ERT survey definition.
%
% Supported types:
%   "single4" : one A-B current pair and one M-N potential pair.
%   "ertLine" : a line/list of electrodes plus quadrupole rows [A B M N].

    arguments
        type string
        options.A (1,3) double = [-5, 0, 0]
        options.B (1,3) double = [ 5, 0, 0]
        options.M (1,3) double = [-2.5, 0, 0]
        options.N (1,3) double = [ 2.5, 0, 0]
        options.Current (1,1) double = 1
        options.Electrodes double = zeros(0, 3)
        options.Quads double = zeros(0, 4)
    end

    type = lower(type);

    switch type
        case "single4"
            electrodes = [options.A; options.B; options.M; options.N];
            quads = [1, 2, 3, 4];

        case "ertline"
            electrodes = options.Electrodes;
            quads = options.Quads;
            if size(electrodes, 2) ~= 3
                error('Electrodes must be an nElectrodes x 3 coordinate array.');
            end
            if size(quads, 2) ~= 4
                error('Quads must be an nMeasurements x 4 array of [A B M N] indices.');
            end

        otherwise
            error('Unsupported survey type: %s', type);
    end

    survey = struct();
    survey.type = type;
    survey.electrodes = electrodes;
    survey.quads = quads;
    survey.current = options.Current;
    survey.nElectrodes = size(electrodes, 1);
    survey.nMeasurements = size(quads, 1);
end
