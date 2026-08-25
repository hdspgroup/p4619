function [sigma, Jsoil, fields] = predict_soil_ec_physics(m, config, aux, ecModel)
%PREDICT_SOIL_EC_PHYSICS Predict EC and local analytic sensitivities.
%
% sigma is in S/m and Jsoil maps stacked property updates to conductivity
% updates at the element level. The inversion vector m lives in transformed
% coordinates for skewed fields, while fields are returned in physical units.

    fields = unpack_soil_fields(m, config);
    uFields = struct();
    for name = config.names
        uFields.(char(name)) = m(config.index.(char(name)));
    end

    allVars = struct();
    for name = config.names
        allVars.(char(name)) = fields.(char(name));
    end
    auxNames = fieldnames(aux);
    for i = 1:numel(auxNames)
        allVars.(auxNames{i}) = aux.(auxNames{i});
    end

    nCells = config.nCells;
    yhat = zeros(nCells, 1);
    dYdU = struct();
    for name = config.names
        dYdU.(char(name)) = zeros(nCells, 1);
    end

    for i = 1:numel(ecModel.termNames)
        term = ecModel.termNames(i);
        coeff = ecModel.coefficients(i);

        if term == "Intercept"
            yhat = yhat + coeff;
            continue;
        end

        if startsWith(term, "state_")
            varName = extractAfter(term, "state_");
            if isfield(allVars, char(varName))
                x = allVars.(char(varName));
                tx = localTransform(varName, x, ecModel);
                dtxdx = localTransformDerivative(varName, x, ecModel);
                dxdu = localInversionDerivative(varName, uFields, config);
                dtdu = dtxdx .* dxdu;
                yhat = yhat + coeff * tx;
                if isfield(dYdU, char(varName))
                    dYdU.(char(varName)) = dYdU.(char(varName)) + coeff * dtdu;
                end
            end
            continue;
        end

        if startsWith(term, "context_")
            varName = extractAfter(term, "context_");
            if isfield(allVars, char(varName))
                x = allVars.(char(varName));
                tx = localTransform(varName, x, ecModel);
                dtxdx = localTransformDerivative(varName, x, ecModel);
                dxdu = localInversionDerivative(varName, uFields, config);
                dtdu = dtxdx .* dxdu;
                yhat = yhat + coeff * tx;
                if isfield(dYdU, char(varName))
                    dYdU.(char(varName)) = dYdU.(char(varName)) + coeff * dtdu;
                end
            end
            continue;
        end

        parts = split(term, "_x_");
        if numel(parts) == 2
            leftName = parts(1);
            rightName = parts(2);
            if isfield(allVars, char(leftName)) && isfield(allVars, char(rightName))
                xl = allVars.(char(leftName));
                xr = allVars.(char(rightName));
                tl = localTransform(leftName, xl, ecModel);
                tr = localTransform(rightName, xr, ecModel);
                dtldx = localTransformDerivative(leftName, xl, ecModel);
                dtrdx = localTransformDerivative(rightName, xr, ecModel);
                dxldu = localInversionDerivative(leftName, uFields, config);
                dxrdu = localInversionDerivative(rightName, uFields, config);
                dtldu = dtldx .* dxldu;
                dtrdu = dtrdx .* dxrdu;
                yhat = yhat + coeff * (tl .* tr);
                if isfield(dYdU, char(leftName))
                    dYdU.(char(leftName)) = dYdU.(char(leftName)) + coeff * (dtldu .* tr);
                end
                if isfield(dYdU, char(rightName))
                    dYdU.(char(rightName)) = dYdU.(char(rightName)) + coeff * (dtrdu .* tl);
                end
            end
        end
    end

    sigma = max(10 .^ yhat, realmin);

    Jblocks = cell(config.nProps, 1);
    for i = 1:config.nProps
        name = config.names(i);
        dSigma = log(10) * sigma .* dYdU.(char(name));
        Jblocks{i} = spdiags(dSigma, 0, nCells, nCells);
    end
    Jsoil = horzcat(Jblocks{:});
end

function tx = localTransform(varName, x, ecModel)
    tx = double(x);
    if any(varName == ecModel.skewedVars)
        tx = log1p(max(tx, 0));
    end
end

function dtx = localTransformDerivative(varName, x, ecModel)
    x = double(x);
    if any(varName == ecModel.skewedVars)
        dtx = 1 ./ (1 + max(x, 0));
        dtx(x < 0) = 0;
    else
        dtx = ones(size(x));
    end
end

function dxdu = localInversionDerivative(varName, uFields, config)
    if isfield(uFields, char(varName))
        dxdu = d_soil_property_d_inversion_var(varName, uFields.(char(varName)), config);
    else
        dxdu = ones(config.nCells, 1);
    end
end
