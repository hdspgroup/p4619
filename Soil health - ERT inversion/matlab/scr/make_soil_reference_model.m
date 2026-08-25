function [mRef, fields] = make_soil_reference_model(config)
%MAKE_SOIL_REFERENCE_MODEL Build a homogeneous starting/reference model.

    fields = struct();
    fields.N = 30 * ones(config.nCells, 1);
    fields.P = 20 * ones(config.nCells, 1);
    fields.K = 120 * ones(config.nCells, 1);
    fields.OC = 2.5 * ones(config.nCells, 1);
    fields.Moisture = 18 * ones(config.nCells, 1);
    fields.Clay = 30 * ones(config.nCells, 1);
    fields.Sand = 45 * ones(config.nCells, 1);
    fields.Silt = 25 * ones(config.nCells, 1);
    mRef = pack_soil_fields(fields, config);
end
