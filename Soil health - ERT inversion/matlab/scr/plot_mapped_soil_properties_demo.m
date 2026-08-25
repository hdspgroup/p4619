function fig = plot_mapped_soil_properties_demo(mesh, fieldsTrue, fieldsMapped, sigmaRecovered)
%PLOT_MAPPED_SOIL_PROPERTIES_DEMO Show staged EC-to-property initialization.

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.04, 0.06, 0.92, 0.82]);
    tl = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Mapped soil properties from recovered EC');

    yMid = round((mesh.gridSize(2) - 1) / 2);

    nexttile; plotField(mesh, sigmaRecovered, yMid, 'Recovered EC (S/m)');
    nexttile; plotField(mesh, fieldsTrue.N, yMid, 'N actual');
    nexttile; plotField(mesh, fieldsMapped.N, yMid, 'N mapped');
    nexttile; plotRatio(mesh, fieldsMapped.N, fieldsTrue.N, yMid, 'N ratio');

    nexttile; plotField(mesh, fieldsTrue.P, yMid, 'P actual');
    nexttile; plotField(mesh, fieldsMapped.P, yMid, 'P mapped');
    nexttile; plotRatio(mesh, fieldsMapped.P, fieldsTrue.P, yMid, 'P ratio');
    nexttile; axis off;

    nexttile; plotField(mesh, fieldsTrue.K, yMid, 'K actual');
    nexttile; plotField(mesh, fieldsMapped.K, yMid, 'K mapped');
    nexttile; plotRatio(mesh, fieldsMapped.K, fieldsTrue.K, yMid, 'K ratio');
    nexttile; plotField(mesh, fieldsTrue.OC, yMid, 'OC actual');

    nexttile; plotField(mesh, fieldsMapped.OC, yMid, 'OC mapped');
    nexttile; plotRatio(mesh, fieldsMapped.OC, fieldsTrue.OC, yMid, 'OC ratio');
    nexttile; plotField(mesh, fieldsMapped.Clay, yMid, 'Clay prior');
    nexttile; plotField(mesh, fieldsMapped.Sand, yMid, 'Sand prior');

    nexttile; plotField(mesh, fieldsMapped.Silt, yMid, 'Silt prior');
    nexttile; plotField(mesh, fieldsMapped.Moisture, yMid, 'Moisture prior');
    nexttile; axis off;
    nexttile; axis off;
end

function plotField(mesh, values, yMid, ttl)
    vals3 = reshape(values, mesh.gridSize - 1);
    imagesc(mesh.x(1:end-1), mesh.z(1:end-1), squeeze(vals3(:, yMid, :))');
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end

function plotRatio(mesh, modelValues, trueValues, yMid, ttl)
    ratio = max(modelValues, realmin) ./ max(trueValues, realmin);
    vals3 = reshape(ratio, mesh.gridSize - 1);
    imagesc(mesh.x(1:end-1), mesh.z(1:end-1), squeeze(vals3(:, yMid, :))');
    set(gca, 'YDir', 'reverse');
    axis tight;
    colorbar;
    title(ttl);
    xlabel('x (m)');
    ylabel('z (m)');
end
