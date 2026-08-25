function fields = evaluate_soil_true_model_at_points(points, refFields, truthSettings)
%EVALUATE_SOIL_TRUE_MODEL_AT_POINTS Evaluate the synthetic truth at xyz points.
%
% points is [n x 3] with columns [x, y, z]. The returned struct contains
% the same physical property fields used by make_soil_true_model.

    if nargin < 3 || isempty(truthSettings)
        truthSettings = localDefaultTruthSettings(refFields);
    else
        truthSettings = localMergeTruthSettings(localDefaultTruthSettings(refFields), truthSettings);
    end

    x = points(:, 1);
    y = points(:, 2);
    z = points(:, 3);

    shallowLayer = exp(-((z - truthSettings.layerCenterZ).^2) / truthSettings.layerWidthDenom);
    anomalyBlob = exp(-((x - truthSettings.anomalyCenterX).^2 + (y - truthSettings.anomalyCenterY).^2) / truthSettings.anomalyXYDenom ...
        - ((z - truthSettings.anomalyCenterZ).^2) / truthSettings.anomalyZDenom);

    fields = refFields;
    baseP = truthSettings.pBackground;
    baseK = truthSettings.kBackground;
    baseOC = truthSettings.ocBackground;
    baseMoisture = truthSettings.moistureBackground;

    fields.N = max(truthSettings.nBackground - truthSettings.nAnomalyAmplitude * anomalyBlob, truthSettings.nAnomalyFloor);
    fields.P = baseP * ones(size(x));
    fields.K = baseK * ones(size(x));
    fields.OC = baseOC * ones(size(x));
    fields.Moisture = baseMoisture * ones(size(x));
    fields.Clay = max(truthSettings.clayBackground + truthSettings.clayLayerAmplitude * shallowLayer, truthSettings.clayFloor);
    fields.Sand = max(truthSettings.sandBackground - truthSettings.sandLayerAmplitude * shallowLayer, truthSettings.sandFloor);
    fields.Silt = 100 - fields.Clay - fields.Sand;
    fields.Silt = max(fields.Silt, truthSettings.siltFloor);
    fields.Sand = 100 - fields.Clay - fields.Silt;
    fields.Sand = max(fields.Sand, truthSettings.sandFloor);
    fields.Silt = 100 - fields.Clay - fields.Sand;
end

function value = localScalarField(x)
    value = x(1);
end

function settings = localDefaultTruthSettings(refFields)
    settings = struct();
    settings.layerCenterZ = 1.9;
    settings.layerWidthDenom = 0.55;
    settings.anomalyCenterX = 0;
    settings.anomalyCenterY = 0;
    settings.anomalyCenterZ = 1.6;
    settings.anomalyXYDenom = 36;
    settings.anomalyZDenom = 0.28;
    settings.nBackground = localScalarField(refFields.N);
    settings.nAnomalyAmplitude = 6;
    settings.nAnomalyFloor = 24;
    settings.pBackground = localScalarField(refFields.P);
    settings.kBackground = localScalarField(refFields.K);
    settings.ocBackground = localScalarField(refFields.OC);
    settings.moistureBackground = localScalarField(refFields.Moisture);
    settings.clayBackground = 20;
    settings.clayLayerAmplitude = 62;
    settings.clayFloor = 5;
    settings.sandBackground = 55;
    settings.sandLayerAmplitude = 36;
    settings.sandFloor = 5;
    settings.siltFloor = 5;
end

function merged = localMergeTruthSettings(defaults, overrides)
    merged = defaults;
    names = fieldnames(overrides);
    for i = 1:numel(names)
        merged.(names{i}) = overrides.(names{i});
    end
end
