function axisVals = make_piecewise_axis(axisLim, segmentEdges, segmentSpacings)
%MAKE_PIECEWISE_AXIS Build a 1D axis with piecewise-uniform spacing.

    arguments
        axisLim (1,2) double
        segmentEdges (1,:) double
        segmentSpacings (1,:) double
    end

    edges = unique([axisLim(1), segmentEdges, axisLim(2)], 'stable');
    if any(diff(edges) <= 0)
        error('Segment edges must be strictly increasing.');
    end
    if numel(segmentSpacings) ~= numel(edges) - 1
        error('Number of segment spacings must equal number of segments.');
    end

    axisVals = edges(1);
    for i = 1:numel(segmentSpacings)
        dx = segmentSpacings(i);
        if dx <= 0
            error('All segment spacings must be positive.');
        end
        a = edges(i);
        b = edges(i + 1);
        seg = a:dx:b;
        if seg(end) ~= b
            seg = [seg, b];
        end
        axisVals = [axisVals, seg(2:end)]; %#ok<AGROW>
    end
    axisVals = unique(axisVals, 'stable');
end
