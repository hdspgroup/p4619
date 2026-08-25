function features = compute_extract_chemistry_features(T, options)
% Compute chemistry features for an aqueous soil-extract EC model.
%
% Expected ion columns are solution concentrations in mg/L:
%   Ca_mgL, Mg_mgL, Na_mgL, K_mgL, Cl_mgL, SO4_mgL, NO3_mgL, HCO3_mgL
%
% Optional soil columns:
%   Clay     clay content in percent
%   OC       organic carbon in g/kg
%   pH_H2O   pH in water
%   pH_CaCl2 pH in CaCl2
%
% If ion chemistry is unavailable, this function falls back to a proxy
% chemistry mode based on pH, extractable K, P, total N, CaCO3, OC, and CEC.
%
% Output fields:
%   sigma_sol_Sm          conductivity estimate or proxy in S/m
%   ionicStrength_molL    measured ionic strength or proxy
%   pH
%   clayFrac
%   OC_gkg
%   CEC
%   CaCO3
%   chemistryMode        1 if ion-based, 0 if proxy-based

    arguments
        T table
        options.TemperatureC double = 25
        options.DefaultPH double = 7
        options.MissingIonAsZero logical = true
        options.ProxyScale_Sm double = 0.1
    end

    ionNames = ["Ca_mgL","Mg_mgL","Na_mgL","K_mgL","Cl_mgL","SO4_mgL","NO3_mgL","HCO3_mgL"];
    molarMass_gmol = [40.078, 24.305, 22.989769, 39.0983, 35.453, 96.06, 62.0049, 61.0168];
    charge = [2, 2, 1, 1, -1, -2, -1, -1];

    % Limiting molar ionic conductivities at 25 C, S cm^2 mol^-1.
    lambda25_Scm2mol = [119.0, 106.1, 50.11, 73.50, 76.35, 160.0, 71.46, 44.5];

    n = height(T);
    conc_mgL = zeros(n, numel(ionNames));
    ionPresent = false(1, numel(ionNames));
    for j = 1:numel(ionNames)
        name = ionNames(j);
        if ismember(name, string(T.Properties.VariableNames))
            conc_mgL(:, j) = T.(name);
            ionPresent(j) = true;
        elseif ~options.MissingIonAsZero
            error('Missing required ion column: %s', name);
        end
    end

    pH = options.DefaultPH .* ones(n, 1);
    if ismember("pH_H2O", string(T.Properties.VariableNames))
        pH = T.pH_H2O;
    elseif ismember("pH", string(T.Properties.VariableNames))
        pH = T.pH;
    elseif ismember("pH_CaCl2", string(T.Properties.VariableNames))
        pH = T.pH_CaCl2;
    end

    clayFrac = zeros(n, 1);
    if ismember("Clay", string(T.Properties.VariableNames))
        clayFrac = T.Clay ./ 100;
    end

    oc = zeros(n, 1);
    if ismember("OC", string(T.Properties.VariableNames))
        oc = T.OC;
    end

    cec = zeros(n, 1);
    if ismember("CEC", string(T.Properties.VariableNames))
        cec = T.CEC;
        cec(~isfinite(cec)) = 0;
    end

    caco3 = zeros(n, 1);
    if ismember("CaCO3", string(T.Properties.VariableNames))
        caco3 = T.CaCO3;
        caco3(~isfinite(caco3)) = 0;
    end

    if any(ionPresent)
        conc_mgL(~isfinite(conc_mgL)) = 0;
        conc_mgL = max(conc_mgL, 0);

        % mg/L is numerically equivalent to g/m^3.
        conc_molm3 = conc_mgL ./ molarMass_gmol;
        conc_molL = conc_molm3 ./ 1000;

        lambda25_Sm2mol = lambda25_Scm2mol .* 1e-4;
        sigma25 = conc_molm3 * lambda25_Sm2mol(:);

        % Common first-order temperature correction for aqueous electrolytes.
        alpha = 0.02;
        sigma_sol = sigma25 .* (1 + alpha .* (options.TemperatureC - 25));
        ionicStrength = 0.5 .* sum(conc_molL .* (charge .^ 2), 2);
        chemistryMode = ones(n, 1);
    else
        % Proxy chemistry mode for bulk soil-property datasets.
        logK = localGetLog1p(T, "K", n);
        logP = localGetLog1p(T, "P", n);
        logN = localGetLog1p(T, "N", n);
        logOC = log1p(max(oc, 0));
        logCEC = log1p(max(cec, 0));
        logCaCO3 = log1p(max(caco3, 0));

        chemIndex = 0.55 .* logK + 0.15 .* logP + 0.10 .* logN + 0.10 .* logCEC + 0.10 .* logCaCO3;
        sigma_sol = options.ProxyScale_Sm .* exp(chemIndex - median(chemIndex, 'omitnan'));
        ionicStrength = 1e-3 .* exp(0.8 .* chemIndex + 0.15 .* (pH - 7) + 0.05 .* logOC);
        chemistryMode = zeros(n, 1);
    end

    features = table(sigma_sol, ionicStrength, pH, clayFrac, oc, cec, caco3, chemistryMode, ...
        'VariableNames', {'sigma_sol_Sm','ionicStrength_molL','pH','clayFrac','OC_gkg','CEC','CaCO3','chemistryMode'});
end

function x = localGetLog1p(T, varName, n)
    x = zeros(n, 1);
    if ismember(varName, string(T.Properties.VariableNames))
        x = T.(varName);
        x(~isfinite(x)) = 0;
        x = log1p(max(x, 0));
    end
end
