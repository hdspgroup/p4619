function electrodes = generate_triangular_planar_electrodes(options)
%GENERATE_TRIANGULAR_PLANAR_ELECTRODES Hexagonal/triangular electrode layout inside a disk.
%
% Produces a regular triangular lattice clipped to a circular footprint.

    arguments
        options.Spacing (1,1) double {mustBePositive} = 5
        options.MaxRadius (1,1) double {mustBePositive} = 25
        options.CenterXY (1,2) double = [0, 0]
        options.Z (1,1) double = 0
        options.RotationDeg (1,1) double = 0
    end

    dx = options.Spacing;
    dy = options.Spacing * sqrt(3) / 2;
    centerXY = options.CenterXY(:).';
    radius = options.MaxRadius;

    rowMax = ceil(radius / max(dy, eps)) + 1;
    xMax = radius + dx;
    pts = zeros(0, 2);

    for row = -rowMax:rowMax
        y = row * dy;
        xOffset = 0.5 * dx * mod(abs(row), 2);
        xVals = (-xMax:dx:xMax) + xOffset;
        yVals = y * ones(size(xVals));
        pts = [pts; [xVals(:), yVals(:)]]; %#ok<AGROW>
    end

    if options.RotationDeg ~= 0
        ang = deg2rad(options.RotationDeg);
        R = [cos(ang), -sin(ang); sin(ang), cos(ang)];
        pts = (R * pts.').';
    end

    radial = hypot(pts(:, 1), pts(:, 2));
    keep = radial <= radius + 1e-9;
    pts = pts(keep, :);

    % Sort by distance to center, then angle, for stable numbering.
    radial = hypot(pts(:, 1), pts(:, 2));
    theta = atan2(pts(:, 2), pts(:, 1));
    [~, order] = sortrows([radial, theta], [1, 2]);
    pts = pts(order, :);

    pts(:, 1) = pts(:, 1) + centerXY(1);
    pts(:, 2) = pts(:, 2) + centerXY(2);
    electrodes = [pts, options.Z * ones(size(pts, 1), 1)];
end
