function model = run_inverse_soil_health_on_data_comp()
% Convenience runner for the inverse model using the current linked dataset.
%
% Inputs:
%   EC, pH, and context variables
% Outputs:
%   N, P, K, CEC, OC

    paths = codex_paths('inverse_soil_health');
    stablePath = fullfile(paths.outputsDir, 'ec_cec_linkage', 'stable_point_table.csv');

    if isfile(stablePath)
        T = readtable(stablePath);
        model = fit_inverse_soil_health_model(T, ...
            'AnalysisName', 'inverse_soil_health', ...
            'ECName', 'EC', ...
            'ConvertECdSmToSm', true);
    else
        s = load(fullfile(paths.dataDir, 'data_comp.mat'));
        varNames = matlab.lang.makeValidName(string(s.varNames), 'ReplacementStyle', 'delete');
        T = array2table(s.data, 'VariableNames', cellstr(varNames));
        model = fit_inverse_soil_health_model(T, ...
            'AnalysisName', 'inverse_soil_health', ...
            'ECName', 'EC', ...
            'ConvertECdSmToSm', true);
    end
end
