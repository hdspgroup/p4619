function m = pack_soil_fields(fields, config)
%PACK_SOIL_FIELDS Pack physical property fields into inversion variables.

    m = zeros(config.totalSize, 1);
    for name = config.names
        m(config.index.(char(name))) = soil_property_to_inversion_var( ...
            name, fields.(char(name)), config);
    end
end
