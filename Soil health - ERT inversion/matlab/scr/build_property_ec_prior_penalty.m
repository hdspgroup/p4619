function [Rec, rhsEc, info] = build_property_ec_prior_penalty(sigmaCurrent, JsoilSigma, prior, weights)
%BUILD_PROPERTY_EC_PRIOR_PENALTY Assemble soft EC-sample prior rows.
%
% This linearizes the point-EC prior directly in conductivity units:
%
%   sigma(m + dm) ~= sigma(m) + Jsigma * dm
%
% using the mapped sample cell for each prior measurement.

    if nargin < 4 || isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'hasMeshMapping') || ~prior.hasMeshMapping
        Rec = sparse(0, size(JsoilSigma, 2));
        rhsEc = zeros(0, 1);
        info = struct('nRows', 0);
        return;
    end

    if ~isfield(prior, 'PointEC') || ~isfield(weights, 'prior') || ~isfield(weights.prior, 'PointEC')
        Rec = sparse(0, size(JsoilSigma, 2));
        rhsEc = zeros(0, 1);
        info = struct('nRows', 0);
        return;
    end

    gamma = weights.prior.PointEC;
    if ~(isnumeric(gamma) && isscalar(gamma) && gamma > 0)
        Rec = sparse(0, size(JsoilSigma, 2));
        rhsEc = zeros(0, 1);
        info = struct('nRows', 0);
        return;
    end

    nRows = prior.nSamples;
    cellIdx = prior.cellIndex(:);
    Rec = sqrt(gamma) * JsoilSigma(cellIdx, :);
    rhsEc = sqrt(gamma) * (prior.PointEC(:) - sigmaCurrent(cellIdx));

    info = struct();
    info.nRows = nRows;
end
