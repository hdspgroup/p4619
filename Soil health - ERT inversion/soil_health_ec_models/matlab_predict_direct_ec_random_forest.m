function [ecHat, details] = predict_direct_ec_random_forest(newData, modelInput)
% Predict EC_sat from a saved direct Random Forest model.
%
% Inputs:
%   newData:
%     - table with variables:
%       pH_CaCl2, pH_H2O, OC, P, N, K, Coarse, Clay, Sand, Silt
%     - or numeric matrix in that exact column order
%   modelInput:
%     - path to MAT file created by train_direct_ec_random_forest
%     - or the loaded model struct itself
%
% Output:
%   ecHat    predicted EC_sat
%   details  struct with transformed/scaled predictors used by the model

    if ischar(modelInput) || isstring(modelInput)
        s = load(modelInput);
        model = s.model;
    else
        model = modelInput;
    end

    predictorNames = model.predictorNames;

    if istable(newData)
        T = newData(:, predictorNames);
    else
        T = array2table(newData, 'VariableNames', cellstr(predictorNames));
    end

    for v = model.skewedPredictors
        T.(v) = log1p(T.(v));
    end

    Xraw = T{:, predictorNames};
    X = (Xraw - model.mu) ./ model.sigma;
    ecHat = exp(predict(model.rf, X)) - 1;

    details = struct();
    details.predictorNames = predictorNames;
    details.transformedPredictors = Xraw;
    details.scaledPredictors = X;
end
