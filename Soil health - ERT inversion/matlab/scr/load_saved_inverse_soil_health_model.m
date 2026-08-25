function model = load_saved_inverse_soil_health_model(modelFile)
%LOAD_SAVED_INVERSE_SOIL_HEALTH_MODEL Load the saved soil-property inverse model.

    s = load(modelFile);
    model = s.model;
end

