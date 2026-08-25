function fig = plot_forward_mesh_structure_3d(mesh, survey)
%PLOT_FORWARD_MESH_STRUCTURE_3D Show the 3D forward mesh geometry for troubleshooting.

    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.10, 0.10, 0.78, 0.72]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Forward mesh structure (%s)', mesh.elementType));

    nodes = mesh.nodes;

    nexttile(1);
    hold on;
    if isfield(mesh, 'meshStyle') && any(strcmpi(string(mesh.meshStyle), ["rect_refined", "octree_style"]))
        localPlotRectMeshShell(mesh);
        scatter3(mesh.elementCenters(:, 1), mesh.elementCenters(:, 2), mesh.elementCenters(:, 3), ...
            10, mesh.elementCenters(:, 3), 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.0);
    elseif isfield(mesh, 'elements') && size(mesh.elements, 2) == 4
        TR = triangulation(mesh.elements, nodes);
        Fb = freeBoundary(TR);
        trisurf(Fb, nodes(:, 1), nodes(:, 2), nodes(:, 3), ...
            'FaceColor', [0.75 0.82 0.92], ...
            'FaceAlpha', 0.18, ...
            'EdgeColor', [0.45 0.55 0.70], ...
            'EdgeAlpha', 0.20);
        scatter3(nodes(:, 1), nodes(:, 2), nodes(:, 3), 8, nodes(:, 3), 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.0);
    else
        scatter3(nodes(:, 1), nodes(:, 2), nodes(:, 3), 8, nodes(:, 3), 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.0);
    end
    if nargin >= 2 && ~isempty(survey) && isfield(survey, 'electrodes')
        scatter3(survey.electrodes(:, 1), survey.electrodes(:, 2), survey.electrodes(:, 3), ...
            36, 'ko', 'LineWidth', 1.2);
    end
    axis equal;
    grid on;
    view(35, 28);
    set(gca, 'ZDir', 'reverse');
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title('3D mesh shell and active cell centers');
    cb = colorbar;
    cb.Label.String = 'z (m)';

    nexttile(2);
    hold on;
    if isfield(mesh, 'meshStyle') && any(strcmpi(string(mesh.meshStyle), ["rect_refined", "octree_style"]))
        plot(mesh.x([1 end end 1 1]), mesh.y([1 1 end end 1]), 'k-', 'LineWidth', 0.8);
        scatter(mesh.elementCenters(:, 1), mesh.elementCenters(:, 2), 18, mesh.elementCenters(:, 3), 'filled');
        zUnique = unique(round(mesh.zCenters, 6));
    else
        zUnique = unique(round(nodes(:, 3), 6));
    end
    zKeep = zUnique(1:min(4, numel(zUnique)));
    if isfield(mesh, 'meshStyle') && any(strcmpi(string(mesh.meshStyle), ["rect_refined", "octree_style"]))
        mask = ismember(round(mesh.elementCenters(:, 3), 6), zKeep);
        scatter(mesh.elementCenters(mask, 1), mesh.elementCenters(mask, 2), 18, mesh.elementCenters(mask, 3), 'filled');
    else
        mask = ismember(round(nodes(:, 3), 6), zKeep);
        scatter(nodes(mask, 1), nodes(mask, 2), 18, nodes(mask, 3), 'filled');
    end
    if nargin >= 2 && ~isempty(survey) && isfield(survey, 'electrodes')
        plot(survey.electrodes(:, 1), survey.electrodes(:, 2), 'ko', 'MarkerSize', 5.5, 'LineWidth', 1.0);
    end
    axis equal;
    grid on;
    xlabel('x (m)');
    ylabel('y (m)');
    title('Plan view of shallow active layers');
    cb2 = colorbar;
    cb2.Label.String = 'z (m)';
end

function localPlotRectMeshShell(mesh)
    x = mesh.x([1 end]);
    y = mesh.y([1 end]);
    z = mesh.z([1 end]);
    [X, Y, Z] = ndgrid(x, y, z);
    corners = [X(:), Y(:), Z(:)];
    edges = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; 4 8; 5 6; 5 7; 6 8; 7 8];
    for i = 1:size(edges, 1)
        pts = corners(edges(i, :), :);
        plot3(pts(:, 1), pts(:, 2), pts(:, 3), '-', 'Color', [0.45 0.55 0.70], 'LineWidth', 0.8);
    end
end
