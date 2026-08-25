clear; clc; close all;

% Demo/control script for the rectilinear FVM ERT forward-modeling pieces:
%   1) refined rectilinear mesh generator
%   2) generic survey struct
%   3) finite-volume operator assembly
%   4) cell-centered forward ERT solve
%   5) diagnostic plotting

saveArg = true;
useParallel = true;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(fullfile(scriptDir, 'scr'));
outDir = fullfile(rootDir, 'outputs', 'ert_forward_demo_fvm_rect');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Setup: mesh

mesh = make_ert_mesh_fvm_rect( ...
    'XLim', [-45, 45], ...
    'YLim', [-45, 45], ...
    'ZLim', [0, 45], ...
    'Spacing', 2.5);

%% Setup: survey

a = 5;
survey = make_ert_survey('single4', ...
    'A', [-a, 0, 0], ...
    'B', [ a, 0, 0], ...
    'M', [-a/2, 0, 0], ...
    'N', [ a/2, 0, 0], ...
    'Current', 1);

%% Setup: resistivity / conductivity model

rho0 = 10;  % ohm m
rhoElements = rho0 * ones(mesh.nElements, 1);

anomaly.center = [2.5, 5, 5];
anomaly.radius = a / 4;
anomaly.rho = rho0 * 100;

d = sqrt(sum((mesh.elementCenters - anomaly.center).^2, 2));
rhoElements(d <= anomaly.radius) = anomaly.rho;
sigmaElements = 1 ./ rhoElements;

%% Assemble and solve

K = assemble_ert_stiffness(mesh, sigmaElements, ...
    'UseParallel', useParallel, ...
    'ElementType', mesh.elementType);

result = solve_ert_forward(mesh, survey, K);

%% Plot

fig = plot_ert_forward_results(mesh, survey, result, ...
    'ReferenceResult', [], ...
    'PropertyElements', rhoElements, ...
    'Anomaly', anomaly, ...
    'SignedLogMinAbs', 1e-3);

if saveArg
    save(fullfile(outDir, 'ert_forward_demo_fvm_rect.mat'), ...
        'mesh', 'survey', 'rhoElements', 'sigmaElements', 'result', 'anomaly');
    exportgraphics(fig, fullfile(outDir, 'ert_forward_demo_fvm_rect.png'), ...
        'Resolution', 220);
end
