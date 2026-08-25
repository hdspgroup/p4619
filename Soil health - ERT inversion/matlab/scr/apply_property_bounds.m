function [mBounded, fieldsBounded, info] = apply_property_bounds(m, config, bounds)
%APPLY_PROPERTY_BOUNDS Clamp the physical property fields to broad bounds.

    fields = unpack_soil_fields(m, config);
    fieldsBounded = fields;
    changed = false;
    nClamped = 0;

    for name = config.names
        key = char(name);
        if ~isfield(bounds, key)
            continue;
        end
        lim = bounds.(key);
        x = fieldsBounded.(key);
        xClamped = min(max(x, lim(1)), lim(2));
        mask = abs(xClamped - x) > 0;
        changed = changed || any(mask(:));
        nClamped = nClamped + nnz(mask);
        fieldsBounded.(key) = xClamped;
    end

    mBounded = pack_soil_fields(fieldsBounded, config);

    info = struct();
    info.changed = changed;
    info.nClamped = nClamped;
end
