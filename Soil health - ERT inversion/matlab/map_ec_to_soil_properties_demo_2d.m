clear; clc; close all;

saveArg = true;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
codexRoot = fileparts(rootDir);
addpath(fullfile(scriptDir, 'scr'));
outDir = fullfile(rootDir, 'outputs', 'map_ec_to_soil_properties_demo_2d');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

ecOnlyFile = fullfile(rootDir, 'outputs', 'ert_ec_only_inversion_demo_2d', 'ert_ec_only_inversion_demo_2d.mat');
if ~isfile(ecOnlyFile)
    error('EC-only 2D inversion result not found: %s', ecOnlyFile);
end

s = load(ecOnlyFile);
inverseModelFile = fullfile(codexRoot, 'Soil health modeling', 'outputs', 'inverse_soil_health', 'inverse_soil_health_model.mat');
inverseModel = load_saved_inverse_soil_health_model(inverseModelFile);

T = table();
T.EC = s.sigmaUpdated(:);
T.pH_CaCl2 = s.aux.pH_CaCl2(:);
T.pH_H2O = s.aux.pH_H2O(:);
T.Clay = s.refFields.Clay(:);
T.Silt = s.refFields.Silt(:);
T.Sand = s.refFields.Sand(:);
T.Coarse = s.aux.Coarse(:);
T.CaCO3 = s.aux.CaCO3(:);

predFields = predict_inverse_soil_health_from_model(T, inverseModel);
fieldsMapped = struct();
fieldsMapped.N = predFields.N;
fieldsMapped.P = predFields.P;
fieldsMapped.K = predFields.K;
fieldsMapped.OC = predFields.OC;
fieldsMapped.Moisture = s.refFields.Moisture;
fieldsMapped.Clay = s.refFields.Clay;
fieldsMapped.Sand = s.refFields.Sand;
fieldsMapped.Silt = s.refFields.Silt;

summary = struct();
summary.targetsEstimated = string(fieldnames(predFields))';
summary.contextCarried = ["Moisture","Clay","Sand","Silt"];

fig = plot_mapped_soil_properties_demo(s.mesh, s.trueFields, fieldsMapped, s.sigmaUpdated);

if saveArg
    save(fullfile(outDir, 'map_ec_to_soil_properties_demo_2d.mat'), ...
        'T', 'fieldsMapped', 'predFields', 'summary', 's');
    exportgraphics(fig, fullfile(outDir, 'map_ec_to_soil_properties_demo_2d.png'), 'Resolution', 220);
end
