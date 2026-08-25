function [ecHat, features] = predict_physics_ec_extract_model(T, modelInput)
% Predict extract EC from a saved physics-informed EC model.
%
% Example:
%   s = load('physics_ec_extract_model.mat');
%   [ecHat, features] = predict_physics_ec_extract_model(newTable, s.model);
%
% This works in either:
%   1) ion-based mode with soluble-ion chemistry, or
%   2) proxy-based mode with bulk soil variables only.

    if ischar(modelInput) || isstring(modelInput)
        s = load(modelInput);
        model = s.model;
    else
        model = modelInput;
    end

    features = compute_extract_chemistry_features(T, ...
        'TemperatureC', model.options.TemperatureC);

    ylog = predict(model.fit, features);
    ecHat = 10 .^ ylog;
end
