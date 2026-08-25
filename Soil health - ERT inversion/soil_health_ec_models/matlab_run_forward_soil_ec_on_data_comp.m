function model = run_forward_soil_ec_on_data_comp()
% Convenience runner for the agreed project architecture:
%   state variables + context variables -> EC
%
% Uses the stable linked point table when available, because that is where
% linked CEC can be combined with EC at the point level.

    paths = codex_paths('forward_soil_ec');
    stablePath = fullfile(paths.outputsDir, 'ec_cec_linkage', 'stable_point_table.csv');

    if isfile(stablePath)
        T = readtable(stablePath);
        model = fit_forward_soil_ec_model(T, ...
            'TargetName', 'EC', ...
            'AnalysisName', 'forward_soil_ec', ...
            'TargetUnits', 'dS/m stored in stable_point_table, converted to S/m', ...
            'ConvertECdSmToSm', true);
    else
        s = load(fullfile(paths.dataDir, 'data_comp.mat'));
        varNames = matlab.lang.makeValidName(string(s.varNames), 'ReplacementStyle', 'delete');
        T = array2table(s.data, 'VariableNames', cellstr(varNames));
        model = fit_forward_soil_ec_model(T, ...
            'TargetName', 'EC', ...
            'AnalysisName', 'forward_soil_ec', ...
            'TargetUnits', 'dS/m from data_comp, converted to S/m', ...
            'ConvertECdSmToSm', true);
    end
end
