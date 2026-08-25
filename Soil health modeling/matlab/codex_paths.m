function paths = codex_paths(analysisName)
% Return standard project paths for the Codex workspace.

    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    [parentDir, thisName] = fileparts(thisDir);

    if strcmpi(thisName, 'matlab')
        projectRoot = parentDir;
        matlabDir = thisDir;
    else
        projectRoot = thisDir;
        matlabDir = fullfile(projectRoot, 'matlab');
    end

    dataDir = fullfile(projectRoot, 'data');
    outputsDir = fullfile(projectRoot, 'outputs');

    if nargin < 1 || strlength(string(analysisName)) == 0
        analysisDir = outputsDir;
    else
        analysisDir = fullfile(outputsDir, char(analysisName));
    end

    dirs = {dataDir, matlabDir, outputsDir, analysisDir};
    for i = 1:numel(dirs)
        if ~exist(dirs{i}, 'dir')
            mkdir(dirs{i});
        end
    end

    paths = struct();
    paths.projectRoot = projectRoot;
    paths.dataDir = dataDir;
    paths.matlabDir = matlabDir;
    paths.outputsDir = outputsDir;
    paths.analysisDir = analysisDir;
end
