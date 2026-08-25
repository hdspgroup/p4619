function ecModel = load_forward_soil_ec_physics_model(coeffFile)
%LOAD_FORWARD_SOIL_EC_PHYSICS_MODEL Load coefficient terms from CSV.

    T = readtable(coeffFile);
    ecModel = struct();
    ecModel.coefficientTable = T;
    ecModel.termNames = string(T.Term);
    ecModel.coefficients = T.MeanCoefficient;
    ecModel.skewedVars = ["N","P","K","CEC","OC","CaCO3","Coarse","Clay","Moisture"];
    ecModel.stateVars = ["N","P","K","OC","Moisture"];
    ecModel.contextVars = ["Clay","Sand","Silt","Coarse","pH_CaCl2","pH_H2O","CaCO3"];
end
