function config = make_soil_property_config(mesh)
%MAKE_SOIL_PROPERTY_CONFIG Define the stacked inversion parameter layout.

    names = ["N","P","K","OC","Moisture","Clay","Sand","Silt"];
    nProps = numel(names);
    nCells = mesh.nElements;

    index = struct();
    for i = 1:nProps
        first = (i - 1) * nCells + 1;
        last = i * nCells;
        index.(char(names(i))) = first:last;
    end

    config = struct();
    config.names = names;
    config.nProps = nProps;
    config.nCells = nCells;
    config.totalSize = nProps * nCells;
    config.index = index;
    config.textureNames = ["Clay","Sand","Silt"];
    config.transformedNames = ["N","P","K","OC","Moisture","Clay","Sand","Silt"];
end
