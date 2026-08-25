% Example workflow for the physics-informed extract EC model.
%
% Replace the variable names in your new data table as needed.
%
% Ion-based version:
%   EC_sat, Ca_mgL, Mg_mgL, Na_mgL, K_mgL, Cl_mgL, SO4_mgL, NO3_mgL,
%   HCO3_mgL, pH_H2O, Clay, OC
%
% Proxy-based version for the data you currently have:
%   EC_sat, Clay, Coarse, Silt, Sand, pH_CaCl2, pH_H2O, OC, CaCO3,
%   N, P, K, and optionally CEC

% T = readtable('your_extract_chemistry_data.csv');
% model = fit_physics_ec_extract_model(T, ...
%     'TargetName', 'EC_sat', ...
%     'TemperatureC', 25, ...
%     'SaveModelFile', 'physics_ec_extract_model.mat');
%
% [ecHat, features] = predict_physics_ec_extract_model(T, model);
% plot(log10(T.EC_sat), log10(ecHat), '.');
% xlabel('log10 measured EC');
% ylabel('log10 predicted EC');
% grid on;
