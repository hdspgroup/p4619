function write_property_inversion_summary(outFile, summary, weights, iterationLog)
%WRITE_PROPERTY_INVERSION_SUMMARY Write a short text summary for the demo.

    if nargin < 4
        iterationLog = struct([]);
    end

    fid = fopen(outFile, 'w');
    fprintf(fid, 'Property-decomposed ERT inversion demo\n');
    fprintf(fid, 'Cells: %d\n', summary.nCells);
    fprintf(fid, 'Measurements: %d\n', summary.nMeasurements);
    fprintf(fid, 'Data misfit before: %.6g\n', summary.dataMisfitBefore);
    fprintf(fid, 'Data misfit after: %.6g\n', summary.dataMisfitAfter);
    fprintf(fid, 'Objective before: %.6g\n', summary.objectiveBefore);
    fprintf(fid, 'Objective after: %.6g\n', summary.objectiveAfter);
    fprintf(fid, 'Update norm: %.6g\n', summary.stepNorm);
    fprintf(fid, 'Last accepted step length: %.6g\n', summary.lastStepLength);
    fprintf(fid, 'Texture-sum residual before: %.6g\n', summary.textureSumResidualBefore);
    fprintf(fid, 'Texture-sum residual after: %.6g\n', summary.textureSumResidualAfter);
    if isfield(summary, 'priorPenaltyAfter')
        fprintf(fid, 'Prior penalty after: %.6g\n', summary.priorPenaltyAfter);
    end
    if isfield(summary, 'priorEcPenaltyAfter')
        fprintf(fid, 'Prior EC penalty after: %.6g\n', summary.priorEcPenaltyAfter);
    end
    if isfield(summary, 'nPriorSamples')
        fprintf(fid, 'Prior samples: %d\n', summary.nPriorSamples);
    end
    if isfield(summary, 'nPriorRowsUsed')
        fprintf(fid, 'Prior rows used: %d\n', summary.nPriorRowsUsed);
    end
    if isfield(summary, 'nPriorEcRowsUsed')
        fprintf(fid, 'Prior EC rows used: %d\n', summary.nPriorEcRowsUsed);
    end
    if isfield(summary, 'priorFile') && ~isempty(summary.priorFile)
        fprintf(fid, 'Prior file: %s\n', summary.priorFile);
    end
    fprintf(fid, 'Iterations completed: %d\n', summary.iterationsCompleted);
    fprintf(fid, 'Exit reason: %s\n', summary.exitReason);
    fprintf(fid, 'Normalized data misfit after: %.6g\n\n', summary.normalizedDataMisfitAfter);
    fprintf(fid, 'Regularization weights\n');
    propNames = fieldnames(weights);
    for i = 1:numel(propNames)
        name = propNames{i};
        if strcmp(name, 'textureSumGamma')
            fprintf(fid, '  textureSumGamma = %.6g\n', weights.textureSumGamma);
        elseif strcmp(name, 'depthBasis')
            fprintf(fid, '  depthBasis.nGaussians = %d\n', weights.depthBasis.nGaussians);
            fprintf(fid, '  depthBasis.includeConstant = %d\n', weights.depthBasis.includeConstant);
        elseif strcmp(name, 'prior')
            priorNames = fieldnames(weights.prior);
            for j = 1:numel(priorNames)
                fprintf(fid, '  prior.%s = %.6g\n', priorNames{j}, weights.prior.(priorNames{j}));
            end
        else
            w = weights.(name);
            fprintf(fid, '  %s: beta=%.6g ax_grad_per_m=%.6g ay_grad_per_m=%.6g az_grad_per_m=%.6g cx_curv_per_m2=%.6g cy_curv_per_m2=%.6g cz_curv_per_m2=%.6g depthBasisGamma=%.6g\n', ...
                name, w.beta, w.ax, w.ay, w.az, w.cx, w.cy, w.cz, w.depthBasisGamma);
        end
    end
    if ~isempty(iterationLog)
        fprintf(fid, '\nIteration log\n');
        for i = 1:numel(iterationLog)
            fprintf(fid, '  iter %d | obj=%.6g | data=%.6g | reg=%.6g | tex=%.6g | prior=%.6g | misfit=%.6g | step=%.6g | alpha=%.3g | accepted=%d\n', ...
                iterationLog(i).iter, iterationLog(i).objective, ...
                iterationLog(i).phiData, iterationLog(i).phiReg, ...
                iterationLog(i).phiTexture, iterationLog(i).phiPrior, iterationLog(i).normalizedDataMisfit, ...
                iterationLog(i).stepNorm, iterationLog(i).stepLength, iterationLog(i).accepted);
        end
    end
    fclose(fid);
end
