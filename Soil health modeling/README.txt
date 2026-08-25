Codex Soil EC Workspace

Project structure

- data
  Stores input datasets such as .mat files.

- matlab
  Stores MATLAB scripts and helper functions.

- outputs
  Stores generated figures, tables, summaries, and saved models.
  Each analysis writes to its own subfolder.

Current data files

- data\data.mat
  Main soil EC dataset used in the earlier EC analyses.

- data\data_comp.mat
  Expanded combined dataset with EC, CEC, repeat-point information, and stable-link proxies.

Main MATLAB scripts

- matlab\soil_ec_analysis.m
  Runs the original soil EC modeling workflow, including plots, model comparison, and sensitivity analysis.
  Outputs go to:
  outputs\soil_ec_analysis

- matlab\direct_ec_model.m
  Compares direct predictive EC models, including polynomial, Random Forest, and boosting approaches.
  Outputs go to:
  outputs\direct_ec

- matlab\train_direct_ec_random_forest.m
  Trains and saves the direct Random Forest EC model.
  Outputs go to:
  outputs\direct_ec_random_forest

- matlab\analyze_ec_cec_linkage.m
  Compares EC prediction with and without linked CEC information using repeat points.
  Outputs go to:
  outputs\ec_cec_linkage

- matlab\fit_physics_ec_extract_model.m
  Fits the physics-informed extract EC model.
  Works in:
  1. ion-based mode if soluble-ion chemistry is available
  2. proxy-based mode if only bulk soil variables are available
  Outputs go to:
  outputs\physics_ec_extract

- matlab\predict_physics_ec_extract_model.m
  Applies a saved physics-informed EC model to new data.

- matlab\compute_extract_chemistry_features.m
  Builds chemistry features for the physics-informed model.

Shared path helper

- matlab\codex_paths.m
  Provides the standard project paths so scripts automatically read from data and write to the correct outputs subfolder.

How to run scripts

Open MATLAB and set the working directory to:

  C:\Users\hdsp\Documents\Data\Codex\matlab

Example commands:

  soil_ec_analysis

  direct_ec_model

  train_direct_ec_random_forest('data.mat', 'direct_ec_rf_model.mat')

  analyze_ec_cec_linkage

Example physics-informed workflow:

  T = readtable('your_data.csv');
  model = fit_physics_ec_extract_model(T, ...
      'TargetName', 'EC_sat', ...
      'TemperatureC', 25, ...
      'SaveModelFile', 'physics_ec_extract_model.mat');
  [ecHat, features] = predict_physics_ec_extract_model(T, model);

Notes

- The scripts expect input .mat files to be in the data folder.
- EC in data_comp.mat is interpreted in dS/m and converted internally to S/m where needed.
- Most model-comparison R^2 values are computed in log10 target space.
- The physics-informed model is currently best treated as a practical template, not a fully mechanistic electrolyte model.
