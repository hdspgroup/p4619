function dxdu = d_soil_property_d_inversion_var(varName, u, config)
%D_SOIL_PROPERTY_D_INVERSION_VAR Derivative of physical x with respect to u.

    u = double(u);
    if is_transformed_soil_property(varName, config)
        dxdu = exp(u);
    else
        dxdu = ones(size(u));
    end
end
