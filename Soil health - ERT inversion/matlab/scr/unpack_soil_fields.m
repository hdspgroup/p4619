function fields = unpack_soil_fields(m, config)
%UNPACK_SOIL_FIELDS Unpack inversion variables into physical property fields.

    fields = struct();
    for name = config.names
        fields.(char(name)) = inversion_var_to_soil_property( ...
            name, m(config.index.(char(name))), config);
    end
end
