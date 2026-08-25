function u = soil_property_to_inversion_var(varName, x, config)
%SOIL_PROPERTY_TO_INVERSION_VAR Map physical property values to inversion vars.

    x = double(x);
    if is_transformed_soil_property(varName, config)
        u = log1p(max(x, 0));
    else
        u = x;
    end
end
