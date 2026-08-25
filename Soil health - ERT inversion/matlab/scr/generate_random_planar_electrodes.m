function electrodes = generate_random_planar_electrodes(options)
%GENERATE_RANDOM_PLANAR_ELECTRODES Random electrode locations inside a disk.

    arguments
        options.NumElectrodes (1,1) double {mustBeInteger, mustBePositive} = 35
        options.MaxRadius (1,1) double {mustBePositive} = 25
        options.CenterXY (1,2) double = [0, 0]
        options.Z (1,1) double = 0
        options.Seed (1,1) double = 20260512
    end

    rng(options.Seed, 'twister');
    theta = 2 * pi * rand(options.NumElectrodes, 1);
    r = options.MaxRadius * sqrt(rand(options.NumElectrodes, 1));

    x = options.CenterXY(1) + r .* cos(theta);
    y = options.CenterXY(2) + r .* sin(theta);
    z = options.Z * ones(options.NumElectrodes, 1);

    electrodes = [x, y, z];
end
