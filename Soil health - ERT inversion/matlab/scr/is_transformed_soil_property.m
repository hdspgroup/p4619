function tf = is_transformed_soil_property(varName, config)
%IS_TRANSFORMED_SOIL_PROPERTY True when the inversion variable is transformed.

    tf = any(string(varName) == config.transformedNames);
end
