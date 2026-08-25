function x = inversion_var_to_soil_property(varName, u, config)
%INVERSION_VAR_TO_SOIL_PROPERTY Map inversion variables back to physical units.

    u = double(u);
    if is_transformed_soil_property(varName, config)
        x = exp(u) - 1;
    else
        x = u;
    end
end
