function [mTrue, fields] = make_soil_true_model(config, mesh, refFields, truthSettings)
%MAKE_SOIL_TRUE_MODEL Create a synthetic heterogeneous soil-property model.

    if nargin < 4
        truthSettings = [];
    end

    fields = evaluate_soil_true_model_at_points(mesh.elementCenters, refFields, truthSettings);

    mTrue = pack_soil_fields(fields, config);
end
