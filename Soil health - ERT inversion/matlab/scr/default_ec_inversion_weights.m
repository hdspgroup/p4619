function weights = default_ec_inversion_weights()
%DEFAULT_EC_INVERSION_WEIGHTS Regularization for EC-only inversion in log10 space.

    weights = struct();
    % For expected layering, favor lateral continuity more than vertical
    % continuity, and use stronger curvature control to suppress
    % checkerboarding without forcing a purely flat model.
    % weights.beta = 3000;
    % weights.ax = 30;
    % weights.ay = 30;
    % weights.az = 30;
    % weights.cx = 200;
    % weights.cy = 200;
    % weights.cz = 60;
    weights.beta = 10^3.5;
    weights.ax = 10^2.5;
    weights.ay = 10^2.5;
    weights.az = 10^2;
    weights.cx = 10^2.5;
    weights.cy = 10^2.5;
    weights.cz = 10^2;

end
