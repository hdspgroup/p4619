function fields = predict_inverse_soil_health_from_model(T, model)
%PREDICT_INVERSE_SOIL_HEALTH_FROM_MODEL Apply saved linear inverse models.

    T = localNormalizeVariableNames(T);
    [X, featureNames] = localInverseFeatures(T, model.sensorVars, model.contextVars);

    fields = struct();
    targetNames = fieldnames(model.targetModels);
    for i = 1:numel(targetNames)
        name = targetNames{i};
        targetModel = model.targetModels.(name);
        coefTbl = targetModel.coefficients;
        intercept = coefTbl.MeanCoefficient(1);
        coefNames = string(coefTbl.Term(2:end));
        coefVals = coefTbl.MeanCoefficient(2:end);
        featureMap = zeros(numel(featureNames), 1);
        for j = 1:numel(featureNames)
            idx = find(coefNames == featureNames(j), 1, 'first');
            if ~isempty(idx)
                featureMap(j) = coefVals(idx);
            end
        end
        yhat = intercept + X * featureMap;
        fields.(name) = 10 .^ yhat;
    end
end

function T = localNormalizeVariableNames(T)
    names = string(T.Properties.VariableNames);
    if ismember("H_2O_Vol_", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2O_Vol_;
    end
    if ismember("H_2OVol", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2OVol;
    end
    if ismember("H_2O_Vol", names) && ~ismember("Moisture", names)
        T.Moisture = T.H_2O_Vol;
    end
end

function [X, names] = localInverseFeatures(T, sensorVars, contextVars)
    X = [];
    names = strings(0, 1);

    for v = sensorVars
        x = localTransform(v, T.(v));
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "sensor_" + v; %#ok<AGROW>
    end

    for v = contextVars
        x = localTransform(v, T.(v));
        X = [X, x]; %#ok<AGROW>
        names(end+1, 1) = "context_" + v; %#ok<AGROW>
    end

    keyContexts = intersect(["Clay","Silt","Sand","CaCO3"], contextVars, 'stable');
    for sv = sensorVars
        xs = localTransform(sv, T.(sv));
        for cv = keyContexts
            xc = localTransform(cv, T.(cv));
            X = [X, xs .* xc]; %#ok<AGROW>
            names(end+1, 1) = sv + "_x_" + cv; %#ok<AGROW>
        end
    end
end

function x = localTransform(varName, x)
    x = double(x);
    skewed = ["EC","EC_sat","N","P","K","CEC","OC","CaCO3","Coarse","Clay","Moisture"];
    if any(varName == skewed)
        x = log1p(max(x, 0));
    end
end

